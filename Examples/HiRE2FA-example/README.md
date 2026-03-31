## Without pip
```bash
### must be run under "HiRE2FA" as current directory
cd ../../HiRE2FA
dir_example=../Examples/HiRE2FA-example/
python3 main.py $dir_example/input/1akx.cg.pdb $dir_example/1akx.fa.pdb
```

## With pip
```bash
### can be run under any current directory
pip install hire2fa
hire2fa input/1akx.cg.pdb 1akx.fa.pdb
```
