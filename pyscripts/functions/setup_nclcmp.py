import matplotlib.colors as mpcrs

def setup_cmap(name,selidx):
    nclcmap='/Users/weiwilliam/AlbanyWork/Utility/colormaps'
    cmapname=name
    f=open(nclcmap+'/'+cmapname+'.rgb','r')
    a=[]
    for line in f.readlines():
        if ('ncolors' in line):
            clnum=int(line.split('=')[1])
        a.append(line)
    f.close()
    b=a[-clnum:]
    c=[]
    for i in selidx[:]:
       if (i==0):
          c.append(tuple(float(y)/255. for y in [255,255,255]))
       elif (i==1):
          c.append(tuple(float(y)/255. for y in [0,0,0]))
       else:
          c.append(tuple(float(y)/255. for y in b[i-2].split()))

    d=mpcrs.LinearSegmentedColormap.from_list(name,c,selidx.size)
    return d
