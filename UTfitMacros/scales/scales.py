import math

varname1 = []
varname2 = []
varlatex = []
int1 = []
int2 = []
int1bis = []
int2bis = []
upperlimit1 = []
upperlimit2 = []

ulowvalue1 = []
ulupvalue1 = []
uloworder1 = []
uluporder1 = []

ulvalue1 = []
ulvalue2 = []
ulorder1 = []
ulorder2 = []
scale1 = []
scale2 = []
n_coeff = 0
n_var = 0
found1 = False
found2 = False

debug = False

data = open('scales.txt')
out = open('scales.tex','w')

coeff = ["ReCK","ImCK","ReCD","ImCD","AbsCBd","AbsCBs"]
latex = [r"$\Re C_K^", r"$\Im C_K^", r"$\Re C_D^", r"$\Im C_D^", r"$|C_{B_d}^", r"$|C_{B_s}^"]
scen = ["GenDF","NMFV"]

indsys = range(0,6)
for inds in indsys :
    indcoeff = range(1,6)
    for indc in indcoeff :
        varname1.append(coeff[inds]+str(indc)+scen[0]+"1")
        varname2.append(coeff[inds]+str(indc)+scen[1]+"1")
        if "Abs" in coeff[inds] :
            varlatex.append(latex[inds]+str(indc)+"|$")
        else :
            varlatex.append(latex[inds]+str(indc)+"$")
            
        if debug :
            print(varname1[n_coeff])
            print(varname2[n_coeff])
            print(varlatex[n_coeff])
        n_coeff = n_coeff+1

if debug :
    print(n_coeff)

loopcoeff = range(0,n_coeff)
for indl in loopcoeff :
    for line in data :
        if "Results for" in line:
            if varname1[indl] in line :
                found1 = True
                if debug :
                    print(line)

                
        if found1 :
            if "at 95.45%" in line:
                found1 = False
                if debug :
                    print(line)
                words = line.split()
                int1.append(words[2])
                if len(words) > 3:
                    int1bis.append(words[3])
                    print("variable "+varname1[indl]+": one splitted interval: you should look into this")

                if debug :
                    print(int1[indl])
                numbers = int1[indl].split(",")
                low = numbers[0].split("[")
                up = numbers[1].split("]")

                lower = float(low[1])
                upper = float(up[0])

                ulow = str(lower).split("e")
                ulowvalue1.append(ulow[0])
                if len(ulow) > 1 :
                    uloworder1.append(ulow[1])
                else :
                    uloworder1.append("0")

                if debug :
                    print(ulow[0])

                ulup = str(upper).split("e")
                ulupvalue1.append(ulup[0])
                if len(ulup) > 1 :
                    uluporder1.append(ulup[1])
                else :
                    uluporder1.append("0")

                if debug :
                    print(ulup[0])

                if abs(lower) > abs(upper):
                    upperlimit1.append(abs(lower))
                    ulvalue1.append(ulow[0])
                    ulorder1.append(ulow[1])            
                else:
                    upperlimit1.append(abs(upper))
                    ulvalue1.append(ulup[0])
                    ulorder1.append(ulup[1])

                if debug :
                    print(ulvalue1[indl])
                    print(ulorder1[indl])

                lam = math.sqrt(1/upperlimit1[indl])
                scale1.append(lam/1000.)
                if debug :
                    print(scale1[indl])

        if "Results for" in line:
            if varname2[indl] in line :
                found2 = True
                if debug :
                    print(line)

                
        if found2 :
            if "at 95.45%" in line:
                found2 = False
                if debug :
                    print(line)
                words = line.split()
                int2.append(words[2])
                if len(words) > 3:
                    int2bis.append(words[3])
                    print("one splitted interval: you should look into this")

                if debug :
                    print(int2[indl])
                numbers = int2[indl].split(",")
                low = numbers[0].split("[")
                up = numbers[1].split("]")

                lower = float(low[1])
                upper = float(up[0])
                if abs(lower) > abs(upper):
                    upperlimit2.append(abs(lower))
                else:
                    upperlimit2.append(abs(upper))

                if debug :
                    print(upperlimit2[indl])

                lam = math.sqrt(1/upperlimit2[indl])
                scale2.append(lam/1000.)
                if debug :
                    print(scale2[indl])

    data.seek(0)


out.write(' Summary \n')
out.write('----------------\n')
indexvar = range(0,n_coeff)
for i in indexvar :
    out.write(varlatex[i])
    out.write('  & $')
    if "Abs" in varname1[i] :
        formatted = "%.1f" %float(ulvalue1[i])
        out.write('< ')
        out.write(str(formatted))
        #out.write(str(ulvalue1[i]))
        out.write(r'\cdot 10^{')
        out.write(str(ulorder1[i]))
        out.write('}$  &  $')
    else :
        if uloworder1[i] != uluporder1[i] :
            diff =  float(uluporder1[i]) - float(uloworder1[i])
            if diff > 0 :
                ulupvalue1[i] = str(float(ulupvalue1[i])*(10*diff))
                uluporder1[i] = uloworder1[i]
            else :
                ulowvalue1[i] = str(float(ulowvalue1[i])*(10*diff))
                uloworder1[i] = uluporder1[i]

        formval1 = "%.1f" %float(ulowvalue1[i])
        formval2 = "%.1f" %float(ulupvalue1[i])
        out.write('[')
        out.write(str(formval1))
        out.write(',')
        out.write(str(formval2))
        out.write('] ')
        #out.write(str(ulvalue1[i]))
        out.write(r'\cdot 10^{')
        out.write(str(uloworder1[i]))
        out.write('}$  &  $')



    if scale1[i] > 1000 :
        formscale1 = scale1[i]/1000
        formatted1 = "%.1f" %float(formscale1)
        out.write(str(formatted1))
        out.write(r'\cdot 10^{3}')        
    else :
        formatted1 = "%.1f" %float(scale1[i])
        out.write(str(formatted1))

    out.write('$  &  $')
    formatted2 = "%.1f" %float(scale2[i])
    out.write(str(formatted2))
    out.write('$ \\\\ \n')

out.write('----------------\n')
out.write(' For the plots \n')
out.write(' General scenario \n')

indsys = range(0,6)
for inds in indsys :
    out.write("Double_t ")
    out.write(coeff[inds])
    out.write("[")
    out.write(str(len(indsys)-1))
    out.write("] = {")
    indcoeff = range(1,6)
    for indc in indcoeff :
        forroot1 = "%.1f" %float(scale1[n_var])
        out.write(str(forroot1))
        if indc != 5 : 
            out.write(", ")
        n_var = n_var+1

    out.write("};\n")

out.write(' NMFV scenario \n')

n_var = 0
indsys = range(0,6)
for inds in indsys :
    out.write("Double_t ")
    out.write(coeff[inds])
    out.write("[")
    out.write(str(len(indsys)-1))
    out.write("] = {")
    indcoeff = range(1,6)
    for indc in indcoeff :
        forroot2 = "%.1f" %float(scale2[n_var])
        out.write(str(forroot2))
        if indc != 5 : 
            out.write(", ")
        n_var = n_var+1

    out.write("};\n")


out.close()


