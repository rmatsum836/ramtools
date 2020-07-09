import numpy as np


o_list = list()
oh_list = list()
f_list = list()
with open('mxene.top', 'r') as f:
    for i, line in enumerate(f):
        if i < 65 or i > 8436:
            continue
        else:
            parse = [s for s in line.split(' ') if s.strip() != '']
            if parse[1] == 'mxene_007':
                index = parse[5]
                o_list.append(int(index))
            elif parse[1] == 'mxene_004':
                index = parse[5]
                oh_list.append(int(index))
            elif parse[1] == 'mxene_006':
                index = parse[5]
                f_list.append(int(index))

init_o = list()
init_oh = list()
init_f = list()
after_o = list()
after_oh = list()
after_f = list()
with open('mxene.gro', 'r') as f:
    for i, line in enumerate(f):
        if i < 3 or i > 402:
            continue
        else:
            parse = [s for s in line.split(' ') if s.strip() != '']
            if int(parse[2]) in o_list:
                zdist = float(parse[-1].split('\n')[0])
                init_o.append(zdist)
            elif int(parse[2]) in oh_list:
                zdist = float(parse[-1].split('\n')[0])
                init_oh.append(zdist)
            elif int(parse[2]) in f_list:
                zdist = float(parse[-1].split('\n')[0])
                init_f.append(zdist)

with open('nvt.gro', 'r') as f:
    for i, line in enumerate(f):
        if i < 3 or i > 402:
            continue
        else:
            parse = [s for s in line.split(' ') if s.strip() != '']
            if int(parse[2]) in o_list:
                zdist = float(parse[5].split('\n')[0])
                after_o.append(zdist)
            elif int(parse[2]) in oh_list:
                zdist = float(parse[5].split('\n')[0])
                after_oh.append(zdist)
            elif int(parse[2]) in f_list:
                zdist = float(parse[5].split('\n')[0])
                after_f.append(zdist)

init_o = np.asarray(init_o)
init_oh = np.asarray(init_oh)
init_f = np.asarray(init_f)
after_o = np.asarray(after_o)
after_oh = np.asarray(after_oh)
after_f = np.asarray(after_f)

dist_o = abs(after_o - init_o)
dist_oh = abs(after_oh - init_oh)
dist_f = abs(after_f - init_f)

print("The mean change of distance for o is: {}\n".format(np.mean(dist_o)))
print("The mean change of distance for oh is: {}\n".format(np.mean(dist_oh)))
print("The mean change of distance for f is: {}".format(np.mean(dist_f)))
