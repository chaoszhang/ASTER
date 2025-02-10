# CASTER Linux walkthrough (Linux版CASTER手把手演示)

[<img src="CASTER.png" width="500"/>](CASTER.png)

# INSTALLATION (安装)

You can choose to clone ASTER repository from GitHub or Gitee: 
(您可以选择从GitHub或Gitee克隆ASTER代码库)

```
git clone https://github.com/chaoszhang/ASTER
```

or (或)

```
git clone https://gitee.com/chaos_zhang/ASTER
```

Now you can enter ASTER directory and install via `make`.
（现在您可以进入ASTER目录并用`make`指令安装）

```
cd ASTER
make
```

Now we enter the directory for our CASTER demo.
(现在我们可以进入CASTER演示用目录了)

```
cd example/caster
```

# Alignments in FASTA format (比对为FASTA格式)

CASTER can take as input a list of alignments in FASTA format. Yes, sequences must be aligned first. Notice that although in this demo I use genes, CASTER actually prefers inter-genic regions.
(CASTER允许将多个FASTA格式比对文件作为输入。没错，你必须先比对好。注意，虽然这个演示用了基因，但CASTER更适合用基因间区域)

```
ls fasta_alignments/*
cat fasta_alignments/gene1.fa
```

## Prepare the input file (准备输入文件)

Put all input FASTA alignment files in a list:
(将所有输入FASTA格式比对文件放到一个列表文件里) 

```
ls fasta_alignments/* > input_fasta_list.txt
cat input_fasta_list.txt
```

Absolute pathes are preferred:
(更推荐用绝对路径)

```
realpath fasta_alignments/* > input_fasta_list.txt
cat input_fasta_list.txt
```

## Run (运行)

Save the output to `caster-site.nw`.
(将输出文件保存到`caster-site.nw`)

```
../../bin/caster-site -i input_fasta_list.txt -o caster-site.nw
cat caster-site.nw
```

Now let's try using 4 threads.
(试着用四个线程看看)

```
../../bin/caster-site -i input_fasta_list.txt -o caster-site.nw -t 4
cat caster-site.nw
```

You can choose a single species as the outgroup.
（您可选择单个物种作为外群）

```
../../bin/caster-site -i input_fasta_list.txt -o caster-site.nw -t 4 --root Orangutan
cat caster-site.nw
```

# Alignments in PHYLIP format (比对为PHYLIP格式)

CASTER can take as input multiple alignments in PHYLIP format. Notice that although in this demo I use genes, CASTER actually prefers inter-genic regions.
(CASTER允许将多个PHYLIP格式比对文件作为输入。注意，虽然这个演示用了基因，但CASTER更适合用基因间区域)

```
ls phylip_alignments/*
cat phylip_alignments/gene1.phy
```

## Prepare the input file (准备输入文件)

Concatenate all input PHYLIP files into one file:
(将所有输入PHYLIP文件合并成一个文件) 

```
cat phylip_alignments/* > concat.phy
cat concat.phy
```

## Run (运行)

Save the output to `caster-site.nw`.
(将输出文件保存到`caster-site.nw`)

```
../../bin/caster-site -i concat.phy -o caster-site.nw
cat caster-site.nw
```

Now let's try using 4 threads.
(试着用四个线程看看)

```
../../bin/caster-site -i concat.phy -o caster-site.nw -t 4
cat caster-site.nw
```

You can choose a single species as the outgroup.
（您可选择单个物种作为外群）

```
../../bin/caster-site -i concat.phy -o caster-site.nw -t 4 --root Orangutan
cat caster-site.nw
```
