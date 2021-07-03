def tc(str):
    return str[0]+str[6:]

def findrk(filename):
    f = open(filename, mode='r')
    data = []
    for line in f:
        line = line.split()
        data.append(line)
    maxrk = 0
    i = 17
    while data[i][0] != 'time':
        rk = int(data[i][2])
        maxrk = max(maxrk, rk)
        i += 1
    return maxrk

def oneD(kernel):
    filename1 = "../comp_1d/Comp_HIFvsHQR_" + kernel + ".txt"
    filename2 = "../comp_1d/Comp_HIFvsHQR_" + kernel + "_4.txt"
    data1 = []
    data2 = []
    f1 = open(filename1, mode='r')
    f2 = open(filename2, mode='r')
    for line in f1:
        line = line.split()
        data1.append(line)
    f1.close()
    for line in f2:
        line = line.split()
        data2.append(line)
    
    ind = [3, 4, 5, 7, 8, 16, 16, 10, 11, 17, 17, 15, 15, 5, 7, 8, 16, 16, 10, 11, 17, 17, 15, 15]
    lines = []
    for i in range(8,20):
        rk1 = findrk("../results_1d/BFF_"+kernel+"_N_"+data1[ind[0]][1]+".txt")
        rk2 = findrk("../results_1d/BFF_"+kernel+"_N_"+data2[ind[0]][1]+"_4.txt")
        lines.append(['$2^{'+str(i)+'}$', tc(data1[ind[1]][3][0:10]), tc(data1[ind[2]][3][0:10]), str(rk1), '1e-3',
                tc(data1[ind[3]][3][0:10]), tc(data1[ind[4]][3][0:10]), data1[ind[5]][9][:-1], tc(data1[ind[6]][11]),
                tc(data1[ind[7]][3][0:10]), tc(data1[ind[8]][3][0:10]), data1[ind[9]][9][:-1], tc(data1[ind[10]][11]),
                data1[ind[11]][10][:-1], tc(data1[ind[12]][12]),
                tc(data2[ind[13]][3][0:10]), str(rk2), '1e-4',
                tc(data2[ind[14]][3][0:10]), tc(data2[ind[15]][3][0:10]), data2[ind[16]][9][:-1], tc(data2[ind[17]][11]),
                tc(data2[ind[18]][3][0:10]), tc(data2[ind[19]][3][0:10]), data2[ind[20]][9][:-1], tc(data2[ind[21]][11]),
                data2[ind[22]][10][:-1], tc(data2[ind[23]][12])])
        for j in range(len(ind)): ind[j] += 18


    for line in lines:
        print('\hline')
        print(line[0]+' & '+line[1]+' & '+line[2]+' & '+line[3]+' & '+line[4]
              +' & '+line[5]+' & '+line[6]+' & '+line[7]+' & '+line[8]+' & '+line[9]
              +' & '+line[10]+' & '+line[11]+' & '+line[12]+' & '+line[13]+' & '+line[14]+'\\\\')
        print('~'+' & '+'~'+' & '+line[15]+' & '+line[16]+' & '+line[17]
              +' & '+line[18]+' & '+line[19]+' & '+line[20]+' & '+line[21]+' & '+line[22]
              +' & '+line[23]+' & '+line[24]+' & '+line[25]+' & '+line[26]+' & '+line[27]+'\\\\')

if __name__ == '__main__':
    oneD('fun_FIO_var4')        