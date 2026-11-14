"""Stream the IEU 1000 Genomes v3 LD reference archive and keep only EUR.bed/.bim/.fam (stops after EUR)."""
import sys, tarfile, urllib.request, os
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
url = "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
want = {"EUR.bed", "EUR.bim", "EUR.fam"}
with urllib.request.urlopen(url) as r:
    tf = tarfile.open(fileobj=r, mode="r|gz")
    for m in tf:
        name = os.path.basename(m.name)
        print("member", m.name, m.size, flush=True)
        if name in want:
            with tf.extractfile(m) as src, open(os.path.join(out, name), "wb") as dst:
                while True:
                    b = src.read(1 << 20)
                    if not b: break
                    dst.write(b)
            want.discard(name)
            print("saved", name, flush=True)
        if not want:
            break
print("done, missing:", want)
