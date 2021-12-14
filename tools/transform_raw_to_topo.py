a = input()
l = [int(x) for x in a.split()]
cnum = len(l)//8
print(cnum)
for i in range(cnum):
    for j in range(8):
        print(l[i*8+j], end='')
        print(' ', end='')
    print('\n', end='')
