import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):

    spec = ([6,16]                               # I6, E16.9
        + [10 for _ in range(24)]               # 24F10.5
        + [13, 13, 13]                          # 3E13.6
        # should the i in the following LC be a _?
        + [12 for _ in range(18)]# for i in (1,12)]#18(1X,E12.5)
        + [9 for _ in range(52)])               # 52F9.5

    df = pd.read_fwf(filename,widths=spec, header=None)
    return df

class Zone:
    def __init__(self, x, y_min, y_max):
        self.x = [x]
        self.y_min = [y_min]
        self.y_max = [y_max]
        self.mid = (y_max + y_min)/2
        
    def append(self, x, y_min, y_max):
        self.x.append(x)
        self.y_min.append(y_min)
        self.y_max.append(y_max)
        self.mid = (y_max + y_min)/2

def calculate_zones(df):

    zones = []
    d = df.iloc[0,9:21]
    conv = d.to_numpy()
    conv = conv[conv >0]
    first_convective = conv[0]
    second_convective = conv[1]

    in_between = d[np.where(d == first_convective)[0][0]+1:np.where(d == second_convective)[0][0]].to_numpy()
    if (in_between < 0).all():
        z = Zone(0, first_convective, second_convective)

    zones.append(z)

    tolerance = 0.1
    q=0
    N = 1
    for count, d in df.iloc[1:,9:21].iterrows():
        first_semi = 0
        tops = []
        bottoms = []
        
        l = d.to_numpy()
        l = np.array(sorted(l, key=lambda x : abs(x)))
        l[l > df.iloc[N, 5]] = 0
        mass = df.iloc[N,5]
        
        conv = l[l > 0]
        
        semi = l[l < 0]
        if len(semi) > 0:
            first_semi = semi[0]
            ind = np.where(first_semi == l)[0][0]
            
            # find bounds around semi-conv
            if ind == 0 or np.sum(l[:ind]) == 0:
                bottoms.append(0)
                tops.append(conv[0])
            else:
                bottoms.append(conv[conv <= np.abs(first_semi)][-1])
                # No top bound above semi present.
                if not (conv > np.abs(first_semi)).any():
                    tops.append(mass)
                else:
                    tops.append(conv[conv > np.abs(first_semi)][0])
                
                # check one below (what if multiple below?)
                if conv[0] != bottoms[-1]:
                    temp =  conv[conv < bottoms[-1]]
                    if len(temp) == 1:
                        bottoms.append(0)
                        tops.append(temp[0])
                    else:
                        bottoms.append(temp[-2])
                        tops.append(temp[-1])

            # look for bounds above the top semi
            if tops[-1] != conv[-1] or tops[-1] != mass:
                temp = conv[conv > tops[-1]]
                if len(temp) == 1:
                    bottoms.append(temp[0])
                    tops.append(mass)
                if len(temp) % 2:
                    for i in range(0,len(temp)-1, 2):
                        bottoms.append(temp[i])
                        tops.append(temp[i+1])
                    bottoms.append(temp[-1])
                    tops.append(mass)
                else:
                    for i in range(0,len(temp), 2):
                        bottoms.append(temp[i])
                        tops.append(temp[i+1])
                        
                        
            for t,b in zip(tops, bottoms):
                if np.abs(t - b) < tolerance:
                    continue
                exists = False
                # check if value within the zone
                for z in zones:
                    if (b <= z.mid and t >= z.mid and np.abs(z.x[-1] - N) < 10):
                        exists = True
                        z.append(N, b, t)
                        break

                # this creates new zones
                if not exists:
                    z = Zone(N, b,t)
                    zones.append(z)

        else:
            # Because we don't have a lot of information to go off, we only expand the old zones
            # there are no semi-convective zones to check what the zone is.
            if len(conv) % 2:
                bottoms.append(0)
                tops.append(conv[0])
                for i in range(1,len(conv)-1,2):
                    bottoms.append(conv[i])
                    tops.append(conv[i+1])
            else:
                for i in range(0,len(conv),2):
                    bottoms.append(conv[i])
                    tops.append(conv[i+1])
            
            
            for t,b in zip(tops, bottoms):
                if np.abs(t - b) < tolerance:
                    continue
                exists = False
                # check if value within the zone
                for z in zones:
                    if (b <= z.mid and t >= z.mid and np.abs(z.x[-1] - N) < 10):
                        exists = True
                        z.append(N, b, t)
                        break
            
        N+=1
    return zones


def get_cores(df):
    return df.iloc[:,5:8].to_numpy().T


def get_index(df):
    return df.iloc[:,0].to_numpy()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
    df = read_file(args.filename)

    zones = calculate_zones(df)
    star = get_cores(df)
    N = get_index(df)
    plt.figure(figsize=(15,15))

    for z in zones:
        plt.fill_between(z.x, z.y_min, z.y_max, color='#d2f8d2')

    plt.plot(N, star[0], lw=5, color='black')
    if (star[1] != 0).any():
        mask = np.where(star[1] != 0)[0][-1]
        star[1][mask:] = star[0][mask:]
    plt.plot(N, star[1], lw=5, color='#1f77b4')
    plt.plot(N, star[2], lw=5, color='#b40424')
    plt.ylim(0,star[0][0]+2)
    plt.xlim(N[0], N[-1])
    plt.savefig('test.png')

