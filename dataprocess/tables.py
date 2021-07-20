def tc(str):
    return str[0:3]+str[6:]

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
    filename1 = "../comp_1d/Comp_HIFvsHQR_" + kernel + "_v4.txt"
    filename2 = "../comp_1d/Comp_HIFvsHQR_" + kernel + "_v6.txt"
    filename3 = "../comp_1d/Comp_HIFvsHQR_" + kernel + "_v8.txt"
    filename4 = "../comp_1d/Comp_HIFvsHQR_" + kernel + "_v12.txt"
    data1 = []
    data2 = []
    data3 = []
    data4 = []
    f1 = open(filename1, mode='r')
    f2 = open(filename2, mode='r')
    f3 = open(filename3, mode='r')
    f4 = open(filename4, mode='r')
    for line in f1:
        line = line.split()
        data1.append(line)
    f1.close()
    for line in f2:
        line = line.split()
        data2.append(line)
    for line in f3:
        line = line.split()
        data3.append(line)
    for line in f4:
        line = line.split()
        data4.append(line)
    
    ind = [3, 4, 5, 7, 8, 16, 16, 10, 11, 17, 17, 15, 15, 5, 7, 8, 16, 16, 10, 11, 17, 17, 15, 15]
    lines = []
    for i in range(8,16):
        rk1 = findrk("../results_1d/BFF_"+kernel+"_N_"+data1[ind[0]][1]+"_v4.txt")
        rk2 = findrk("../results_1d/BFF_"+kernel+"_N_"+data2[ind[0]][1]+"_v6.txt")
        rk3 = findrk("../results_1d/BFF_"+kernel+"_N_"+data3[ind[0]][1]+"_v8.txt")
        lines.append(['$2^{'+str(i)+'}$', tc(data1[ind[1]][3][0:10]), tc(data1[ind[2]][3][0:10]), str(rk1), '1e-4',
                tc(data1[ind[3]][3][0:10]), tc(data1[ind[4]][3][0:10]), data1[ind[5]][9][:-1], tc(data1[ind[6]][11]),
                tc(data1[ind[7]][3][0:10]), tc(data1[ind[8]][3][0:10]), data1[ind[9]][9][:-1], tc(data1[ind[10]][11]),
                data1[ind[11]][10][:-1], tc(data1[ind[12]][12]),
                tc(data2[ind[13]][3][0:10]), str(rk2), '1e-6',
                tc(data2[ind[14]][3][0:10]), tc(data2[ind[15]][3][0:10]), data2[ind[16]][9][:-1], tc(data2[ind[17]][11]),
                tc(data2[ind[18]][3][0:10]), tc(data2[ind[19]][3][0:10]), data2[ind[20]][9][:-1], tc(data2[ind[21]][11]),
                data2[ind[22]][10][:-1], tc(data2[ind[23]][12]),
                tc(data3[ind[13]][3][0:10]), str(rk3), '1e-8',
                tc(data3[ind[14]][3][0:10]), tc(data3[ind[15]][3][0:10]), data3[ind[16]][9][:-1], tc(data3[ind[17]][11]),
                tc(data3[ind[18]][3][0:10]), tc(data3[ind[19]][3][0:10]), data3[ind[20]][9][:-1], tc(data3[ind[21]][11]),
                data3[ind[22]][10][:-1], tc(data3[ind[23]][12])])
                
        if i < 15:
            rk4 = findrk("../results_1d/BFF_"+kernel+"_N_"+data4[ind[0]][1]+"_v12.txt")
            line = [tc(data4[ind[13]][3][0:10]), str(rk4), '1e-12',
                tc(data4[ind[14]][3][0:10]), tc(data4[ind[15]][3][0:10]), data4[ind[16]][9][:-1], tc(data4[ind[17]][11]),
                tc(data4[ind[18]][3][0:10]), tc(data4[ind[19]][3][0:10]), data4[ind[20]][9][:-1], tc(data4[ind[21]][11]),
                data4[ind[22]][10][:-1], tc(data4[ind[23]][12])]
            lines[-1] = lines[-1] + line

        for j in range(len(ind)): ind[j] += 18


    for i in range(8,16):
        line = lines[i-8]
        print('\hline')
        print(line[0]+' & '+line[1]+' & '+line[2]+' & '+line[3]+' & '+line[4]
              +' & '+line[5]+' & '+line[6]+' & '+line[7]+' & '+line[8]+' & '+line[9]
              +' & '+line[10]+' & '+line[11]+' & '+line[12]+' & '+line[13]+' & '+line[14]+'\\\\')
        print('~'+' & '+'~'+' & '+line[15]+' & '+line[16]+' & '+line[17]
              +' & '+line[18]+' & '+line[19]+' & '+line[20]+' & '+line[21]+' & '+line[22]
              +' & '+line[23]+' & '+line[24]+' & '+line[25]+' & '+line[26]+' & '+line[27]+'\\\\')
        print('~'+' & '+'~'+' & '+line[28]+' & '+line[29]+' & '+line[30]
              +' & '+line[31]+' & '+line[32]+' & '+line[33]+' & '+line[34]+' & '+line[35]
              +' & '+line[36]+' & '+line[37]+' & '+line[38]+' & '+line[39]+' & '+line[40]+'\\\\')
        if i < 15:
            print('~'+' & '+'~'+' & '+line[41]+' & '+line[42]+' & '+line[43]
              +' & '+line[44]+' & '+line[45]+' & '+line[46]+' & '+line[47]+' & '+line[48]
              +' & '+line[49]+' & '+line[50]+' & '+line[51]+' & '+line[52]+' & '+line[53]+'\\\\')

def cond():
    filename = "../condnum.txt"
    data = []
    f = open(filename, mode='r')
    for line in f:
        line = line.split()
        data.append(line)
    
    lines = []
    ind1 = [76, 22, 40, 58, 4]
    ind2 = [77, 23, 41, 59, 5]
    for i in range(8, 13):
        lines.append(['$2^{'+str(i)+'}$', '$A$', data[ind1[0]][4], data[ind1[1]][4], data[ind1[2]][4], data[ind1[3]][4], data[ind1[4]][4], '$A^{*}A$', data[ind2[0]][4], data[ind2[1]][4], data[ind2[2]][4], data[ind2[3]][4], data[ind2[4]][4]])
        for j in range(len(ind1)):  ind1[j] += 3
        for j in range(len(ind2)):  ind2[j] += 3

    for line in lines:
        print('\hline')
        print(line[0]+' & '+line[1]+' & '+line[2]+' & '+line[3]+' & '+line[4]+' & '+line[5]+' & '+line[6]+'\\\\')
        print('~'+' & '+line[7]+' & '+line[8]+' & '+line[9]+' & '+line[10]+' & '+line[11]+' & '+line[12]+'\\\\')


def sb1d(kernel):
    filename1 = "../comp_1d/Regularization_" + kernel + "_L1.txt"
    filename2 = "../comp_1d/Regularization_" + kernel + "_TV-L1.txt"
    filename3 = "../comp_1d/Regularization_" + kernel + "_GS-TV-L1.txt"
    data1 = []
    data2 = []
    data3 = []
    f1 = open(filename1, mode='r')
    f2 = open(filename2, mode='r')
    f3 = open(filename3, mode='r')
    for line in f1:
        line = line.split()
        data1.append(line)
    for line in f2:
        line = line.split()
        data2.append(line)
    for line in f3:
        line = line.split()
        data3.append(line)
    ind = [3, 4, 5, 7, 8, 14]
    lines = []
    for i in range(8,17):
        lines.append(['$2^{'+str(i)+'}$', tc(data1[ind[1]][3][0:10]), tc(data1[ind[2]][3][0:10]), '1e-7',
                tc(data1[ind[3]][3][0:10]), tc(data1[ind[4]][3][0:10]), data1[ind[5]][11][0:2-1], tc(data1[ind[5]][13]),
                data2[ind[5]][11][0:-1], tc(data2[ind[5]][13]), data3[ind[5]][11][0:-1], tc(data3[ind[5]][13])])
        for i in range(len(ind)): ind[i] += 16 
    
    for i in range(8,17):
        line = lines[i-8]
        print('\hline')
        print(line[0]+' & '+line[1]+' & '+line[2]+' & '+line[3]+' & '+line[4]
              +' & '+line[5]+' & '+line[6]+' & '+line[7]+' & '+line[8]+' & '+line[9]+' & '+line[10]+' & '+line[11]+'\\\\')


if __name__ == '__main__':
    # oneD('fun_FIO_5')      
    # cond()
    sb1d('fun_FIO_5')