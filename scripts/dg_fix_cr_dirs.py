import os,sys
root=sys.argv[1]; nf=nd=0
for d in sorted(os.listdir(root)):
    if "\r" not in d: continue
    clean=d.replace("\r","")
    if os.path.exists(os.path.join(root,clean)): sys.exit("collision "+clean)
    p=os.path.join(root,d)
    for f in os.listdir(p):
        if "\r" in f:
            os.rename(os.path.join(p,f), os.path.join(p,f.replace("\r",""))); nf+=1
    os.rename(p, os.path.join(root,clean)); nd+=1
print("renamed dirs",nd,"files",nf)
