import math, random
import numpy as np, pandas as pd

#run CDC's import data R script to get transformed_data.csv for import into this program
#this program will run the same heuristic algorithm for generating a distance matrix as euk_heuristic R scrip
#but uses pairwise deletion instead of imputing missing values

hapfile = "transformed_data_66B.csv"
outfile = "heuristic_distmat_66B.csv"

NULLVAL = "NA"
locinames = ("AA","AC","AD","AE","AF","AG","AH","AJ","AK","AL","AM","AO_PART_A","AO_PART_B","AP","AQ","AR","AS","AU","AV","AW","AX","AY","AZ","CA","CB","CC","CD","X360i2_PART_A","X360i2_PART_B","X378_PART_A","X378_PART_B","MSR_PART_A","MSR_PART_B","FA","FB","FC","FD","FE","FF","FG","FH","FI","FL","FM","FP","FQ","FR","FT","FU","FV","FW","MTJ")
## provide ploidy of each locus - ordered the same as the locinames variable
ploidyval = (2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1)
ploidy = { locinames[i]:ploidyval[i] for i in range(len(ploidyval))}

data = pd.read_csv(hapfile)
ids = data['ids'].tolist()
#df = data.drop('ids', 1)
locicolumns = list(data.columns)
nids = len(ids)
nloci = len(locinames)
print(ids)

frequencies = dict() #dictionary of dictionaries
observeddatamatrix = dict() #dictionary of lists
H_nu = dict()
for loc in locinames:
    observeddatamatrix[loc] = []
    sub = data.filter(regex = loc)
    #sub = sub.fillna('X')
    sub2 = sub.to_numpy()
    freq1 = dict()
    for c in range(len(sub2[0])):
        count = 0
        l1 = []
        for r in range(len(sub2)):
            item = sub2[r][c]
            if type(item) == str: #remove nan
                freq1[item] = freq1.get(item, 0) + 1
                l1.append(item)
                count += 1
            else:
                l1.append(NULLVAL)
        if count > 0:
            observeddatamatrix[loc].append(l1)
    total = sum(freq1.values(), 0.0)
    h_sum = 0.0
    for k,v in freq1.items():
        freq1[k] = v / total
        h_sum -= freq1[k] * math.log2(freq1[k])
    frequencies[loc] = freq1
    H_nu[loc] = h_sum

m = {k:1 for k in locinames}
def sub_per_locus(isolate1,isolate2,loc):
    v1 = [ col[isolate1] for col in observeddatamatrix[loc] if col[isolate1] != NULLVAL]
    v2 = [ col[isolate2] for col in observeddatamatrix[loc] if col[isolate2] != NULLVAL] 
    if len(v1) == 0 or len(v2) == 0:
        return NULLVAL
    p1 = [ frequencies[loc][x] for x in v1]
    p2 = [ frequencies[loc][x] for x in v2]

    def fnz(z, jj):
        total = 0
        for ii in range(jj, jj*2+1):
            total -= ii*(z == ii)
        return total

    s1 = set(v1)
    s2 = set(v2)
    shared_alleles = s1.intersection(s2)
    x = len(s1) + len(s2)
    y = len(shared_alleles)
    if len(shared_alleles) > 0:
        temp_shared = [0 for _ in range(nids)]
        notmissing = [0 for _ in range(nids)]
        for i in range(nids):
            all_alle = set()
            for col in observeddatamatrix[loc]:
                if col[i] != NULLVAL:
                    notmissing[i] += 1
                    all_alle.add(col[i])
            for alle in shared_alleles:
                if alle in all_alle:
                    temp_shared[i] += 1
        #print (temp_shared)
        #print (notmissing)
        P = 0.0
        denom = 0.0
        for i in range(nids):
            if temp_shared[i] == len(shared_alleles):
                P += 1.0
            if notmissing[i] > 0:
                denom += 1.0
        P = P * P / denom / denom      
    else:
        P = 1.0
    k = 1 * (y == 0) + P * (y > 0)
    if ploidy[loc] > 1:
        n = min(len(s1),len(s2))
        w = x * (n > 1) + 4 * (n == 1) * (x == 2) + (1+x) * (n== 1) *  (x > 2)
        jj = 2 * (m[loc] == 1) + m[loc] * (m[loc] > 1)
        z = 3 * (((2*(n == 1) + 1 * (y == 1) * (x > 2)))==3) + 2 * (y *(n>1)*(jj>=y)+jj*(n > 1)* (y > jj) + jj * (n==1) * (x==2) * (y==1))
        delta_nu_raw = w * (y == 0) + 2 * jj * (y > 0) + fnz(z, jj)
        delta_nu = H_nu[loc] * ( delta_nu_raw * (delta_nu_raw > 0) + P * (delta_nu_raw == 0))*k
        delta = delta_nu
    else:
        delta_ex_raw = 2 * x * (y == 0)
        delta_ex = H_nu[loc] * ( delta_ex_raw * (delta_ex_raw > 0) + P * (delta_ex_raw == 0))*k
        delta = delta_ex
    return delta

#sub_per_locus(0,4,"AA")

def pairwisedistance_heuristic(isolate1,isolate2):
    sumdelta = 0.0
    num_used = 0.0
    ok = True
    for loc in locinames:
        d = sub_per_locus(isolate1,isolate2,loc)
        #delta.append(d)
        try:
            sumdelta += d
            num_used += 1.0
        except TypeError:
            ok = False
    return sumdelta * nloci / num_used

matrix = [[0 for _ in range(nids)] for _ in range(nids)]
for i in range(nids):
    for j in range(i + 1, nids):
        matrix[i][j] = matrix[j][i] = pairwisedistance_heuristic(i,j)

#output 
out_file = open(outfile, 'w')
out = "Seq_ID," + ",".join(ids) + "\n"
out_file.write(out)
for i in range(nids):
    out = ids[i] + "," + ",".join(map(str,matrix[i])) + "\n"
    out_file.write(out)
    #print(c[1], c[0])
out_file.close()



