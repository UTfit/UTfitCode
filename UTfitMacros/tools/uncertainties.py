import math

var = []
mean = []
err = []
relerr = []

debug = False

data = open('measurements.conf')
out = open('uncertainties.txt','w')

coeff = ["Vcb", "Vub", "Epsk", "Dmd", "Dms", "S2b", "S2a", "Gamma", "BRbtaunuexp"]

loopcoeff = range(0,len(coeff))
for indl in loopcoeff :
    for line in data :
        if not "//measurement" in line:
            if "measurement" in line:
                #print(line)
                if coeff[indl] in line :
                    words = line.split()
                    #print(words)
                    var.append(words[1])
                    val1 = words[3].split(",")
                    mean.append(val1[0])
                    val2 = words[4].split(",")
                    err.append(val2[0])
                    if float(val1[0]) != 0 :
                        relerr.append(float(val2[0])*100./float(val1[0]))
                    else :
                        relerr.append(0.)
    data.seek(0)

out.write(' Summary \n')
out.write('----------------\n')
indexvar = range(0,len(var))
for i in indexvar :
    out.write(str(i))
    out.write('/ ')
    out.write(var[i])
    out.write(': ')
    out.write(mean[i])
    out.write(', ')
    out.write(err[i])
    out.write(' -> ')
    formatted = "%.1f" %relerr[i]
    out.write(str(formatted))
    out.write('%\n')


out.write('----------------\n')

ind1 = 6
ind2 = 2
out.write(var[ind1]+"/"+var[ind2])
out.write(': ')
val1 = float(mean[ind1])
err1 = float(err[ind1])
val2 = float(mean[ind2])
err2 = float(err[ind2])
ratio = val1/val2
out.write(str(ratio))
out.write(', ')
ratioerr = ratio*math.sqrt(pow(err1/val1,2)+pow(err2/val2,2));
out.write(str(ratioerr))
out.write(' -> ')
ratiorelerr = ratioerr*100/ratio
formatted = "%.1f" %ratiorelerr
out.write(str(formatted))
out.write('%\n')

ind1 = 10
ind2 = 11
out.write(var[ind1]+"/"+var[ind2])
out.write(': ')
val1 = float(mean[ind1])
err1 = float(err[ind1])
val2 = float(mean[ind2])
err2 = float(err[ind2])
ratio = val1/val2
out.write(str(ratio))
out.write(', ')
ratioerr = ratio*math.sqrt(pow(err1/val1,2)+pow(err2/val2,2));
out.write(str(ratioerr))
out.write(' -> ')
ratiorelerr = ratioerr*100/ratio
formatted = "%.1f" %ratiorelerr
out.write(str(formatted))
out.write('%\n')

out.write('----------------\n')
out.close()


