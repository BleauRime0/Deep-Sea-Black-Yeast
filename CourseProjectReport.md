# 深海黑酵母*Hortaea werneckii* MC743的基因组组装、注释和比较分析

王帜飞，牛潘豪

[TOC]

###  摘要

***Hortaea werneckii*** 是一种极端耐盐的黑色酵母，属于灰孢霉目，最近在地中海的不同站点和深度中被分离出来，并证明是优势真菌物种。为了探索这些地中海分离株的基因组特征，我们小组在阅读前人的研究后决定沿用文献中提及的数据，选择合适的生物信息工具软件，在 Linux 系统下以文献中 SRA 数据库的数据为基础，尝试对海洋***H. werneckii*** 菌株进行基因组组装并进一步进行比较基因组学分析，最终与文献得到的结果进行对比分析，总结在复刻实验中存在的不足与需要改进之处。

关键词：***Hortaea werneckii***；黑酵母；基因组组装；基因注释；比较基因组学

### 前言

 ***Hortaea werneckii*** 属于黑色酵母，这是一个多源的多态和黑色真菌群体，在许多情况下呈现出多极端耐受的生活方式。迄今为止，该物种主要是因其显著的耐盐性而受到研究，它是唯一已知能够在广泛的盐度范围内生长的真菌，从0 M NaCl到饱和的5.1 M NaCl。因此，该真菌被用作研究真核细胞在高盐度条件下的渗透适应策略和耐受性所涉及的分子机制的模式生物。

***H. werneckii***是一种全球适应性广泛的真菌物种，适应于世界各地的高盐度环境，如富营养太阳盐场中的盐卤，这被认为是其主要的生态位。事实上，已经证明在NaCl浓度高于20％的这些栖息地中，它是优势真菌物种，因此其存在已报道于不同的咸水环境，如海水、海滩土壤、盐场微生物垫、浸泡木材、鱼类、珊瑚和盐沼植物。迄今为止，已经对来自不同来源和不同国家的12个 ***H. werneckii*** 菌株的基因组进行了测序。

本研究的目的是对来自地中海 Vector 和 Geostar 海洋站的两个 ***H. werneckii*** 菌株的完整基因组进行全新的基因组组装和注释，并进行比较基因组学分析，分别命名为 **MC848** 和 **MC873**。



## 数据集与方法

​        根据前言所示，本项目的工作分为七个步骤：原始测序质量评估，基因组组装，重复序列和可移动元件预测，基因和RNA预测，基因注释和翻译，比较基因组学分析和蛋白质功能注释。

### 一. 原始测序质量评估

####         *本板块用 bash脚本对测序数据量进行统计，并创建了工作环境使用 FastQC工具对测序文件进行分析，从而生成质控报告。*

在序列读段归档 (Sequence Reads Archive, SRA) 数据库中，获取生物项目识别符**PRJNA641248**的原始测序数据[Run Selector :: NCBI (nih.gov)](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA641248&o=acc_s%3Aa)。下载 BioSample [SAMN15347964](https://www.ncbi.nlm.nih.gov/biosample/SAMN15347964) 包含的 SRR12072898, SRR12072899, SRR12072900 共3个Run。

#### 1. 测序数据量统计

编写bash脚本，对于目录下的测序原始文件`*.fastq`，分析其 Total reads 和 Read length，结果输出到命令行。使用附录中的`bash`脚本`fastq_analysis.sh`统计测序数据。

由于全部3个测序文件较大，此处仅以`SRR12072898.fastq`的输出作为示例：

```
Analyzing file: SRR12072899.fastq
Total Reads: 61670952tq
Read Length Range: 150 - 150
```

统计结果如下：**Read Length**: 150 bp

| Sample      | SRR12072898 | SRR12072899 | SRR12072900 |
| ----------- | ----------- | ----------- | ----------- |
| Total Reads | 15,417,738  | 27,512,404  | 30,177,094  |

#### 2. 原始数据碱基质量分布

1. 创建虚拟环境`bio2503`，安装**FastQC**。

   ```bash
   conda create -n bio2503
   conda activate bio2503
   sudo apt update
   sudo apt install fastqc
   ```

2. 将测序文件`SRR12072898.fastq.gz`放在`quality_control`目录下，使用FastQC对测序文件进行分析，并生成质控报告。

   ```bash
   cd ~/quality_control
   fastqc SRR12072898.fastq.gz
   ```

3. FastQC生成`SRR12072898_fastqc.html`和`SRR12072898_fastqc.zip`，分别包含质控结果的详细信息。在浏览器中打开HTML文件查看结果。重复以上步骤对其余2个测序文件进行质检，测序质量均符合要求。

   ![](C:\Users\HataYou\Pictures\Saved Pictures\1717476760645.jpg)



### 二. 基因组组装流程

####       *本板块使用 `Trimmomatic`工具对测序数据进行修剪，并选择 `SPAdes`工具作为基因组装的工具，最终使用 QUAST与 Pilon实现了组装后的处理。*

#### 1. 数据预处理

1. **质量控制**：使用`FastQC`检查测序数据的质量。

   创建指定的输出目录`~/qc_reports/`：

   ```bash
   mkdir -p ~/qc_reports/
   ```

    使用`FastQC`检查测序数据的质量：

   ```bash
   fastqc SRR12072898.fastq -o qc_reports/
   ```

2. **除低质量数据和接头**：使用`Trimmomatic`进行数据的修剪。

   使用 `conda` 安装`Trimmomatic`:

   ```bash
   conda install -c bioconda trimmomatic
   ```

   前往 `Trimmomatic` 的 GitHub 页面下载最新的版本：

   ```bash
   wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
   ```

   解压下载的文件：

   ```bash
   unzip Trimmomatic-0.39.zip
   ```

   在安装和配置 `Trimmomatic` 后，运行 `Trimmomatic`进行数据修剪：

   ```bash
   trimmomatic SE -phred33 SRR12072898.fastq SRR12072898_trimmed.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
   ```

#### 2. 基因组组装

1. **选择适合的组装工具，比如`SPAdes`、`ABySS`、`SOAPdenovo`等。这里选择的是`SPAdes`**：

   安装`SPAdes`：

   ```bash
   conda install -c bioconda spades
   ```

   使用`SPAdes`进行基因组组装：

   ```bash
   spades.py --s1 SRR12072898_trimmed.fastq -o spades_output --isolate -m 16 -t 4
   ```

   其中，`-m 16` 表示使用最多 16GB 内存，`-t 4` 表示使用 4 个线程。原因是`SPAdes` 运行时遇到了内存不足的问题，需要通过设置 `SPAdes` 使用较低的 `k-mer` 或调整参数来降低内存消耗。

   运行界面如下图所示：

   ![image-20240604132618575](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240604132618575.png)

#### 3. 组装后处理

   1. **组装质量评估**：使用QUAST评估组装的质量。

      使用`conda`安装 QUAST：

      ```bash
      conda install -c bioconda quast
      ```

      在安装的过程中，遇到了多个依赖包的问题，导致无法安装 QUAST。为解决问题，选择重新创建新的`conda`环境并安装了必要的依赖：

      ```bash
      conda create -n quast_env python=3.7
      conda activate quast_env
      conda install -c bioconda -c conda-forge numpy matplotlib
      ```

      从 QUAST GitHub 页面下载并安装 QUAST：

      ```bash
      wget https://github.com/ablab/quast/releases/download/quast_5.0.2/quast-5.0.2.tar.gz
      tar -xzf quast-5.0.2.tar.gz
      cd quast-5.0.2
      ./install_full.sh
      ```

      将 QUAST 添加到系统路径中：

      ```bash
      export PATH=$PATH:/home/niupanhao/bio2502/project/quast-5.0.2
      ```

      运行 QUAST:

      ```bash
      quast.py spades_output/contigs.fasta -o quast_output
      ```

      运行途中，遇到Python 3.8 及更高版本中 `cgi.escape` 被移除的问题，最终通过修改QUAST源代码来修复：

      ```bash
      vim quast_libs/site_packages/jsontemplate/jsontemplate.py
      cgi.escape
      #替换为
      html.escape
      #在文件上方添加
      import html
      ```

      运行得到结果显示了 contig 数量和长度分布，以及基因组的全长：

      ![image-20240605192310151](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240605192310151.png)

2. **基因组抛光**：

   生成BAM文件：

   ```bash
   # 以BWA为例
   bwa index spades_output/contigs.fasta
   bwa mem spades_output/contigs.fasta reads_1.fastq reads_2.fastq > aligned.sam
   # 将SAM文件转换为BAM文件并进行排序
   samtools view -bS aligned.sam > aligned.bam
   samtools sort aligned.bam -o aligned_sorted.bam
   # 索引BAM文件
   samtools index aligned_sorted.bam
   ```

   使用 Pilon 进行基因组抛光: 

      ```bash
   java -jar pilon-1.24.jar --genome spades_output/contigs.fasta --frags aligned_sorted.bam --output polished_assembly
      ```

   基因组组装最终得到的结果为`contigs.fasta`文件，截取部分为如下图。将基因组组装文件命名为`MC743.fna`，用于下游分析。

   ![image-20240604132545661](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240604132545661.png)




### 三. 重复序列和可移动元件预测

####         *本版块安装 `RepeatMasker`和 `RepeatModeler`并配置相应环境，进行重复序列数据库的构建，并列出各种类型的可移动元件和重复序列的数量、长度和占比。*

#### 1. RepeatMasker

1. 创建`conda `环境`repeatmakser`，进入环境中安装`python`、`perl`和相关的`python`包。

   ```bash
   conda create -n repeatmasker
   conda activate repeatmasker
   conda install python=3.7
   conda install perl
   conda install h5py
   ```

2. 安装预编译的`rmblastn`，并测试是否可用。

   ```bash
   wget -c https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz
   tar -zxvf rmblast-2.13.0+-x64-linux.tar.gz
   rmblastn -h
   ```

3. 安装**TRF**，软件已编译好，下载后修改权限。

   ```bash
   wget -c https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
   ln -s ~/trf409.linux64 /home/data/user222/biosoft/trf
   chmod 755 trf
   ```

4. 下载**RepeatMasker**，然后进入目录配制数据集。

   ```bash
   wget -c http://repeatmasker.org/RepeatMasker/RepeatMasker-4.1.4.tar.gz
   tar -zxvf RepeatMasker-4.1.4.tar.gz
   cd RepeatMasker
   ```

   下载**Dfam**数据库：

   ```bash
   cd ~/RepeatMasker/Libraries
   wget https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz
   wget https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz.md5sum
   gunzip Dfam.h5.gz
   ```

   下载**RepBase**，解压后自动添加到已有的`Libraries`中：

   ```bash
   wget https://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20181026.tar.gz
   tar -zxvf RepBaseRepeatMaskerEdition-20181026.tar.gz
   ```

5. 编译`Dfam`和`repbase`中的重复序列，并配置工具。

   ```bash
   perl ./configure 
   ```

   将软件软连接到环境变量，先查看当前的环境变量：

   ```bash
   which python
   ln -s ~/RepeatMasker/RepeatMasker ~/anaconda3/envs/repeatmasker/bin/python
   ```

#### 2. RepeatModeler

1. 安装**recon**并编译：

   ```bash
   http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
   tar -zxvf RECON-1.08.tar.gz
   cd RECON-1.08/src
   make
   make install
   ```

2. 安装**RepeatScout**并编译：

   ```bash
   wget -c http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
   tar -zxvf RepeatScout-1.0.6.tar.gz
   cd RepeatScout-1.0.6
   make
   ```

3. 安装需要的**perl**包：

   ```bash
   conda install -c bioconda perl-json
   conda install -c bioconda perl-file-which
   conda install -c bioconda perl-uri
   conda install -c bioconda perl-devel-size
   conda install -c eumetsat perl-lwp
   ```

4. 安装**cd-hit**并编译：

   ```bash
   wget -c https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
   tar -zxvf cd-hit-v4.8.1-2019-0228.tar.gz
   cd cd-hit-v4.8.1-2019-0228
   make
   cd cd-hit-auxtools
   make
   ```

5. 安装**mafft**，编译并修改`Makefile`：

   ```bash
   wget -c https://mafft.cbrc.jp/alignment/software/mafft-7.505-without-extensions-src.tgz
   tar -zxvf mafft-7.505-without-extensions-src.tgz
   cd mafft-7.505-without-extensions
   cd core
   ```

   修改`Makefile`中的`PREFIX`如下：

   ```makefile
   PREFIX = /usr/local
   # 改为
   PREFIX = ~/mafft-7.505-without-extensions
   ```

   退出，重新编译：

   ```bash
   make clean
   make
   make install
   ```

6. 安装下列**UCSC** tools，无需编译，放在目录`~/biolab/UCSCTOOLS`下：

   ```bash
   wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
   wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
   wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
   ```

7. 安装**RepeatModeler**，配置依赖路径：

   ```bash
   wget -c http://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.4.tar.gz
   tar -zxvf RepeatModeler-2.0.4.tar.gz
   cd RepeatModeler-2.0.4
   perl ./configure # 根据提示依次填写各依赖项路径
   ```

8. 将 RepeatModeler 添加至环境变量：

   ```bash
   vim ~/.bashrc
   # 在文件末尾添加
   export PATH="$PATH:~/RepeatModeler"
   ```

#### 3. 重复序列注释

1. 创建目录`~/biolab/transposons`，使用`MC743.fna`基因组，为RepeatModeler创建索引数据。

   ```bash
   BuildDatabase -name mc743 -engine ncbi MC743.fna
   ```

   `-engine ncbi`表示使用 rmblast，`-name mc743`表示数据库的名字为 mc743.

2. 运行 RepeatModeler 从头预测：

   ```bash
   RepeatModeler -database mc743 -engine ncbi -threads 20 &> mc743.out &
   ```

   `-database` 和上一步一致，`-engine` 和上一步一致，`-threads` 表示线程数。

3. 运行时间相对比较久，过程中的临时文件存放在`RM4111468_.XXX`文件夹下。

   运行结束后，得到了`mc743-families.fa`和`ath-families.stk`。 前者是找到的重复序列，后者是Stockholm格式的种子联配文件(seed alignment file)。此外，`MC743.fna.tbl`直观列出了各种类型的可移动元件和重复序列的数量、长度和占比。

<img src="C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191612419.png" alt="image-20240606191612419" style="zoom: 50%;" />

### 四. 基因/rRNA/tRNA预测

####         *本版块分别针对基因 /rRNA/tRNA，通过不同的工具 `GeneMark`/ `RNAmmer`/ `tRNAscan-SE `进行基因预测*。

#### 1. 真核基因预测

1. **软件安装**：GeneMark官网 ([GeneMark™ download (gatech.edu)](http://topaz.gatech.edu/GeneMark/license_download.cgi)) 勾选 GeneMark-ES/ET/EP + ver 4.72_lic，LINUX 64 kernel 3.10 - 5，下载程序和密钥：

   >  Please download program [here](http://topaz.gatech.edu/GeneMark/tmp/GMtool_H9xVq/gmes_linux_64_4.tar.gz)
   >
   >  Please download key [32_bit](http://topaz.gatech.edu/GeneMark/tmp/GMtool_H9xVq/gm_key_32.gz) or [64_bit](http://topaz.gatech.edu/GeneMark/tmp/GMtool_H9xVq/gm_key_64.gz)

2. **环境配置**：把`gm_key_64`和`gmes_linux_64_4`的压缩包放在家目录下解压。

   编辑环境变量：

   ```bash
   sudo vim ~/.bashrc
   ```

   在文末添加以下内容：

   ```bash
   export PATH=~/biolab/gmes_linux_64_4:$PATH
   export GM_KEY=~/biolab/gm_key_64
   ```

   保存并退出，然后重新加载配置文件：

   ```bash
   source ~/.bashrc
   ```


3. **基因预测**：进入 Genemark 程序目录`~/gmes_linux_64_4`，创建新的工作目录，防止运行时冲突，然后进入工作目录：

   ```bash
   mkdir ~/gmes_linux_64_4/work
   cd ~/gmes_linux_64_4/work
   ```

   把密钥移动到家目录下作为不可见文件：

   ```bash
   mv ~/gm_key_64 ~/.gm_key_64
   ```

   对基因组文件`MC743.fna` （路径为`~/gene_predict/MC743.fna`）进行预测：

   ```bash
   perl ~/gmes_linux_64_4/gmes_petap.pl --sequence ~/gene_predict/MC743.fna --ES --cores 4
   ```


> Perl 相关问题：执行最后一步时可能缺少依赖模块，例如`Hash::Merge`、`MCE::Mutex`、`YAML`等。`YAML`用`sudo apt-get install libyaml-perl`安装，其余用`cpan`内置`install`安装。多次运行`perl`指令以检查模块有无缺漏。

<img src="C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191507840.png" alt="image-20240606191507840" style="zoom: 50%;" />

#### 2. rRNA预测

1. 在**RNAmmer**官网（[RNAmmer 1.2 - DTU Health Tech - Bioinformatic Services](https://services.healthtech.dtu.dk/services/RNAmmer-1.2/)）填写信息，申请试用权限，下载软件包`rnammer-1.2.Unix.tar.gz`。移动到目录`~/biolab/`下并解压：

   ```bash
   tar -xzvf rnammer-1.2.Unix.tar.gz
   ```

2. 在**HMMER**官网（[HMMER](http://www.hmmer.org/download.html)）下载最新版本软件包，移动到相同目录下并解压。

   ```bash
   tar -xzvf hmmer-3.4.tar.gz
   ```

   解压后在hmmer-2.3.2文件夹下安装流程如下：

   ```bash
   cd hmmer-3.4
   ./configure
   make
   make check
   make install
   ```

3. **环境配置**：进入RNAmmer的目录并配置环境变量。确保系统中已经安装了依赖的 Perl 模块。

   ```bash
   vim rnammer
   ```

   设置安装途径和`HMMSEARCH`途径：

   ```ini
   my $INSTALL_PATH = "~/rnammer-1.2";
   $HMMSEARCH_BINARY = "~/hmmer-3.4/src/hmmsearch";
   ```

   将`rnammer`脚本的执行权限设置为可执行：

   ```bash
   chmod +x rnammer
   ```

4. **基因预测**：设置适当的参数，运行 RNAmmer 进行rRNA预测：

   ```bash
   ./rnammer -S euk -m lsu,ssu,tsu -gff MC743_rRNA.gff < ../gene_predict/MC743.fna
   ```

   - `-S bac` 指定物种类型为真核生物（`euk`），可根据样本调整为`arc`（古菌）或`bak`（细菌）。

   - `-m lsu,ssu,tsu` 指定要预测的rRNA类型，分别代表大亚基、小亚基和转运RNA。
   - `-gff MC743_rRNA.gff` 指定输出文件格式为GFF，并将结果保存到`../gene_predict/MC743_rRNA.gff`文件中。`< ../gene_predict/MC743.fna` 指定输入基因组文件路径。

> Perl 相关问题：执行最后一步时可能缺少依赖模块，例如`XML::Simple`可用`cpan`内置`install`安装。

![image-20240606191436528](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191436528.png)

#### 3. tRNA预测

1. 从**tRNAscan-SE**官网（[tRNAscan-SE Search Server (ucsc.edu)](http://lowelab.ucsc.edu/tRNAscan-SE/)）下载软件包`trnascan-se-2.0.12.tar.gz`。移动到目录`~/biolab/`下并解压：

   ```bash
   tar -zvxf trnascan-se-2.0.12.tar.gz
   ```

2. 进入解压后的目录并安装`tRNAscan-SE`，依次运行配置和安装步骤：

   ```bash
   cd tRNAscan-SE-2.0.12
   ./configure
   make
   make install
   ```

3. **环境配置**：确保`tRNAscan-SE`的可执行文件在`PATH`中。

   ```bash
   export PATH=~/tRNAscan-SE-2.0.12/bin:$PATH
   ```

   将上述命令添加到`~/.bashrc`，以便登录时自动生效：

   ```bash
   echo 'export PATH=~/tRNAscan-SE-2.0.12/bin:$PATH' >> ~/.bashrc
   source ~/.bashrc
   ```

4. **安装依赖项**：从**Infernal**官网（[Infernal: inference of RNA alignments (eddylab.org)](http://eddylab.org/infernal/)）下载

   ```bash
   cd ~/biolab
   wget http://eddylab.org/infernal/infernal-1.1.5.tar.gz
   tar -xzvf infernal-1.1.5.tar.gz
   cd infernal-1.1.5
   ./configure --prefix=$HOME/local
   make
   make install
   ```

   将`Infernal`的安装路径添加到`PATH`和`~/.bashrc`中。

   ```bash
   export PATH=$HOME/local/bin:$PATH
   echo 'export PATH=$HOME/local/bin:$PATH' >> ~/.bashrc
   source ~/.bashrc
   ```

   `which cmsearch`和`which cmscan`确认已正确安装并在`PATH`中。输出`cmsearch`的可执行文件路径为`~/local/bin/cmsearch`，由于`tRNAscan-SE`的默认路径与此不相同，创建符号链接，将实际的`cmsearch`路径链接到`tRNAscan-SE`寻找的默认路径。

   ```bash
   sudo ln -s ~/local/bin/cmsearch /usr/local/bin/cmsearch
   sudo ln -s ~/local/bin/cmscan /usr/local/bin/cmscan
   ```

5. **基因预测**：设置适当的输入和输出路径，运行`tRNAscan-SE`进行tRNA预测：

   ```bash
   cd ~/gene_predict
   tRNAscan-SE -o MC743_tRNAscan.txt MC743.fna
   ```

   `-o MC743_tRNAscan_results.txt` 指定输出文件名，将结果保存到`MC743_tRNAscan_results.txt`中。

经过以上三步，得到基因组`MC743.fna`的基因预测文件`MC743.genemark.gtf`，rRNA注释文件`MC743_rRNA.gff`和tRNA注释文件`MC743_tRNAscan.txt`。

![image-20240606191356894](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191356894.png)



### 五. 基因注释和翻译

####         *本板块使用 `gffread`工具得到 CDSs核苷酸序列，并使用 `transeq` 工具将其翻译成蛋白质序列。*

#### 1. 蛋白质编码序列 (CDSs)

使用**GeneMark-ES**对一个基因组组装文件`MC743.fna`进行了基因预测，得到了基因组注释文件`MC743.genemark.gtf`，可以给出基因组里蛋白质编码序列（CDSs）在基因组某个scaffold上的位置，并标明正义链或反义链，但是并不能直接给出CDSs的核苷酸序列。`gffread`可以从GFF/GTF文件中提取序列，是`Cufflinks`软件包的一部分，也可以单独使用。

1. 激活虚拟环境`bio2503`，在`~/biolab/`目录下创建工作目录`gene_annotation`并将`MC743.gtf`和`MC743.fna`复制到该目录下。

2. **软件安装**：使用`apt-get`安装`gffread`

   ```bash
   sudo apt install gffread
   ```

   使用`conda`安装`samtools`：

   ```bash
   conda install -c bioconda samtools
   ```

3. **文件预处理**：由于`gffread`工具不接受注释行有空格的情况，对`GTF`和`fna`文件的注释部分进行处理，保留`scaffold`编号作为识别特征。使用附录中的`python`脚本`extract_scaffold.py`提取特征。

4. **解析基因组**：使用`gffread`提取CDSs序列：

   ```bash
   gffread MC743.gtf -g MC743.fna -x MC743.cds.fasta
   ```

   `MC743.genemark.gtf` 是GTF注释文件，`MC743.fna` 是您的基因组文件。`-x MC743_CDS.fasta` 指定输出文件名，将CDSs保存到`MC743.cds.fasta`中。

   经过以上步骤，得到了编码序列文件`MC743.cds.fasta`，包含18954条核苷酸序列。

<img src="C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191200441.png" alt="image-20240606191200441" style="zoom: 50%;" />

#### 2. 核苷酸序列翻译

由于编码基因在COG和KEGG数据集上的注释，通常是基于氨基酸序列而非核苷酸序列的**blastp**比对，我们需要把获得的CDSs翻译成多肽序列，以便进行下一步的功能分析。在 Linux 操作系统上，可以使用生物信息学软件工具如 `EMBOSS` 将核苷酸序列翻译成蛋白质序列。

1. 安装 `EMBOSS` 工具包：

   ```bash
   sudo apt-get install emboss
   ```

2. **翻译**核苷酸序列：CDS 核苷酸序列文件位于 `~/biolab/gene_annotation/MC743.cds.fa`。使用 `transeq` 工具将其翻译成蛋白质序列。

   ```bash
   transeq -sequence MC743.cds.fasta -outseq MC743.protein.fa
   ```

经过以上步骤，得到了从**MC743**基因组预测出的18954条多肽序列文件`MC743.protein.fa`。

<img src="C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191103757.png" alt="image-20240606191103757" style="zoom:50%;" />



### 六. 比较基因组学分析

####         *本版块使用 OrthoMCL软件进行氨基酸序列的比对与相似性聚类，并使用 MUMmer软件进行两个基因组序列的共线性分析。*

#### 1. 同源基因分析

采用**OrthoMCL**软件对所有参与分析的物种的氨基酸序列进行比对，选取一定阈值（在30%～80%之间，视具体情况而定）进行相似性聚类，获得同源基因的列表。统计每一个蛋白聚类cluster的物种分布情况，可用于种内的泛基因组、核心基因组的研究。

1. **软件依赖**：

   - BLAST来进行序列比对，通过以下命令安装BLAST+

   ```bash
   sudo apt-get update
   sudo apt-get install ncbi-blast+
   ```

   - MCL（Markov Cluster Algorithm）用于聚类分析

   ```bash
   sudo apt-get install mcl
   ```

   - OrthoMCL需要MySQL数据库来存储和管理数据

   ```bash
   sudo apt-get install mysql-server
   ```

2. 下载和安装**OrthoMCL**

   ```bash
   wget https://orthomcl.org/common/downloads/software/v2.0/orthomclSoftware-v2.0.9.tar.gz
   tar -xzvf orthomclSoftware-v2.0.9.tar.gz
   cd orthomclSoftware-v2.0.9
   export PATH=$PATH:~/orthomclSoftware-v2.0.9/bin
   ```

3. 创建一个用于OrthoMCL的数据库和用户：

   ```bash
   sudo mysql -u root -p
   ```

   在MySQL提示符下输入以下命令：

   ```sql
   CREATE DATABASE orthomcl;
   CREATE USER 'orthomcl'@'localhost' IDENTIFIED BY 'password';
   GRANT ALL PRIVILEGES ON orthomcl.* TO 'orthomcl'@'localhost';
   FLUSH PRIVILEGES;
   EXIT;
   ```

4. 在`orthomclSoftware-v2.0.9`目录下，复制示例配置文件并进行编辑：

   ```bash
   cp doc/orthomcl.config.template orthomcl.config
   nano orthomcl.config
   ```

   编辑`orthomcl.config`文件，更改部分内容如下：

   ```ini
   dbVendor=mysql
   dbConnectString=dbi:mysql:orthomcl:localhost:mysql_socket=/var/run/mysqld/mysqld.sock
   dbLogin=orthomcl
   dbPassword=password
   orthomclInstallSchema
   ```

   安装数据库模式：

   ```bash
   orthomclInstallSchema orthomcl.config
   ```

5. 将FASTA文件格式化为OrthoMCL所需的格式：

   ```bash
   mkdir -p ~/homology/formatted
   orthomclAdjustFasta MC743 ~/homology/MC743.protein.fa 1
   orthomclAdjustFasta MC848 ~/homology/MC848.protein.fa 1
   mv *.fasta ~/homology/formatted/
   ```

6. 生成BLAST数据库和执行BLAST

   ```bash
   makeblastdb -in ~/homology/formatted/MC743.fasta -dbtype prot
   makeblastdb -in ~/homology/formatted/MC848.fasta -dbtype prot
   
   blastp -query ~/homology/formatted/MC848.fasta -db ~/homology/formatted/MC743.fasta -outfmt 6 -evalue 1e-5 -out ~/homology/blast.out
   ```

   加载 blast 结果

   ```bash
   orthomclBlastParser ~/homology/blast.out ~/homology/formatted > ~/homology/similarSequences.txt
   ```

7. 过滤相似序列：编辑过滤配置文件

   ```bash
   cp doc/orthomclFilterFasta.template orthomclFilterFasta.config
   vim orthomclFilterFasta.config
   ```

   修改配置文件中的`similarSequences`路径

   ```ini
   similarSequencesFile=~/homology/similarSequences.txt
   ```

   然后执行过滤

   ```bash
   orthomclLoadBlast orthomcl.config ~/homology/similarSequences.txt
   orthomclPairs orthomcl.config orthomclPairs.log cleanup=yes
   orthomclDumpPairsFiles orthomcl.config
   ```

8. 对生成的`mclInput`文件，使用MCL进行聚类分析。

   ```bash
   mcl mclInput --abc -I 1.5 -o mclOutput
   ```

   这里的`-I 1.5`是膨胀参数，根据实际情况调整来优化聚类结果。

9. 将MCL的输出转换为OrthoMCL的群集文件格式。

   ```bash
   orthomclMclToGroups OG 1 < mclOutput > homoGeneCluster.txt
   ```

   `OG`是群集名称的前缀，`1`是起始编号。`homoGeneCluster.tsv`文件包含了最终的同源基因群集信息。

   ![image-20240608170705029](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240608170705029.png)

#### 2. 共线性分析

采用**MUMmer**软件，进行两个基因组序列的共线性分析，可以从宏观清晰地显示基因组水平上的插入、缺失、翻转、易位等现象。

1. **软件安装**：先使用以下命令安装依赖包：

   ```bash
   sudo apt-get install gcc make zlib1g-dev libncurses-dev
   ```

   前往 MUMmer 官方页面（[Release Mummer 4.0.0 release candidate 1 · mummer4/mummer (github.com)](https://github.com/mummer4/mummer/releases/tag/v4.0.0rc1)）下载 MUMmer 4.0.0rc1 的源代码压缩包并解压。

   ```bash
   tar -xzf mummer-4.0.0rc1.tar.gz
   ```

   进入 MUMmer 源代码目录，配置和编译 MUMmer，然后安装：

   ```bash
   cd mummer-4.0.0rc1
   ./configure
   make
   sudo make install
   ```

2. 设置**环境变量**，确保 MUMmer 的二进制文件可以在终端中被找到。

   ```bash
   export PATH=~/mummer-4.0.0rc1:$PATH
   ```

3. 创建工作目录`~/biolab/collinearity/`，将两个基因组文件复制到该目录下，使用 `nucmer` 对两个基因组序列进行比对，`-p MC743_vs_MC848` 指定输出文件的前缀：

   ```bash
   nucmer -p MC743_vs_MC848 MC743.fna MC848.fna
   ```

4. **过滤数据**：为提高比对结果的有效性和可视化图像的美观度，保留相似度大于 80% 和长度大于 1000 bp 的比对：

   ```bash
   delta-filter -1 -i 90 -l 1000 MC743_vs_MC848.delta > filtered.delta
   ```

   `-i`表示相似度 (identity)，`-l`表示长度 (length)，`-1`只保留一对一的比对（每个查询序列只比对到一个目标序列）。

5. 使用 `show-coords` 工具生成比对结果的坐标文件：

   ```bash
   show-coords -rcl filtered.delta > MC743_vs_MC848.coords
   ```

   `-r` 按参考序列排序，`-c`显示覆盖率，`-l` 输出更详细的比对信息。

6. **可视化**比对结果：MUMmer 提供了多种可视化工具，可使用 `mummerplot` 生成可视化图：

   ```bash
   mummerplot --png -p MC743_vs_MC848 filtered.delta
   ```

   `-p MC743_vs_MC848` 指定输出文件的前缀。进入生成的`MC743_vs_MC848.gp`文件修改前缀：

   ```bash
   vim MC743_vs_MC848.gp
   # 前缀修改为
   set terminal postscript
   set output 'MC743_vs_MC848.ps'
   # 保存并退出vim
   gnuplot MC743_vs_MC848.gp
   ```

   使用`ghostscript`命令行工具，将`MC743_vs_MC848.ps`转换成`collinearity.pdf`文件。

   ```bash
   gs -sDEVICE=pdfwrite -o collinearity.pdf MC743_vs_MC848.ps
   ```

   <img src="C:\Users\HataYou\Pictures\Saved Pictures\collinearity_00.png" style="zoom: 50%;" />



### 七、蛋白质功能注释

#### *本版块在本地建立 COG数据集，使用 BLAST和 HMMER进行蛋白质功能注释，并比较 MC743和 MC848基因功能的差异。*

#### COG功能分析

**COG**是Clusters of Orthologous Groups of proteins的缩写 (http://www.ncbi.nlm.nih.gov/COG/)。COG是在对已完成基因组测序的物种的蛋白质序列进行相互比较的基础上构建的，COG数据库选取的物种包括各个主要的系统进化谱系。每个COG家族至少由来自3个系统进化谱系的物种的蛋白所组成，所以一个COG对应于一个古老的保守结构域。构成每个COG的蛋白被假定来自于同一个祖先蛋白。进行COG数据库比对可以对预测蛋白进行功能注释、归类以及蛋白进化分析。

在Linux操作系统上进行蛋白质功能注释可以通过使用COG数据库和一些生物信息工具来实现。

1. **软件安装**：使用`conda`安装**BLAST**和**HMMER**：

   ```bash
   conda install -c bioconda blast
   conda install -c bioconda hmmer
   ```

2. **COG**数据库：访问[Index of /pub/COG/COG2020/data (nih.gov)](https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/)以获取最新的COG数据库文件，下载并解压：

   ```bash
   wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/COG2020.fa.gz
   gunzip cog-20.fa.gz
   ```

3. 构建**BLAST数据库**，将蛋白质序列与COG数据库进行比对：

   ```bash
   makeblastdb -in cog-20.fa -dbtype prot -out cog-20
   blastp -query MC743.protein.fa -db cog-20 -out MC743_COGs_blast.txt -evalue 1e-5 -outfmt 6
   ```

   亦可使用**隐马尔可夫模型**比对：

   - 下载并解压PFAM数据库（包含隐马尔可夫模型）

   ```bash
   wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz
   gunzip Pfam-A.hmm.gz
   ```

   - 构建HMMER数据库

   ```bash
   hmmpress Pfam-A.hmm
   ```

   - 使用 HMMER 进行比对

   ```bash
   hmmscan --tblout MC743_COGs_hmmer.txt Pfam-A.hmm MC743.protein.fa
   ```

4. 经过对比，HMMER算法在本地的效率远高于BLAST。根据HMMER的比对结果，附录中的`python`脚本`hmmer_results.py`用于过滤低相似度比对，解析并注释蛋白质功能，输出表格`protein_function.tsv`。

   以上方法过滤 **E value** > 1e-10 的比对，保留基因名称(target_name)、基因编号(accession)、蛋白质序列编号(query_name)、E 值(evalue)、得分(score)、偏差(bias)和蛋白质功能(description)。MC743 的18954个蛋白质序列共获得18792条高质量注释，MC848 的18962个蛋白质序列共获得18865条高质量注释，两者基因数量和功能广度未体现明显差异（相对差异 < 0.1%）。下面是 `MC743_COGs_hmmer.txt` 和 `MC743_protein_function.tsv` 的截取部分。

![image-20240606191725075](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191725075.png)

![image-20240606191937022](C:\Users\HataYou\AppData\Roaming\Typora\typora-user-images\image-20240606191937022.png)



### 讨论

在本项目中，我们根据文献所给途径下载了一株来自3400米深度（MC873，正文内容中误记为MC743，下文沿用此编号）的菌株的测序文件，使用 Linux 上的生信工具重新组装，并沿用作者提供的来自2500米深度(MC848)的菌株现成组装，进行基因组注释和基因功能等分析。在原文中，研究者获得了MC848和MC873基因组的高质量组装(>99.8%)，并证明它们是二倍体。对于*H. werneckii* 的两个基因组，检测到了高水平的杂合性(约4.6%)。这是二倍体物种和/或经历过全基因组复制(WGD)的真核生物谱系的典型特征，由两个单倍体祖先之间的偶发的种内杂交事件导致。该研究中进行的系统发生学分析显示存在两个独立的聚类群。有趣的是，从地中海不同站点的3400米和2500米深处回收的菌株聚集在一起，表明这些基因型存在一定的环境特异性。

文章以列表的形式给出了原始测序文件的读长分布和质量得分分布，但是在 NCBI 上发布的 SRA 版本中已经将读长统一修剪为150 bp，质量统一控制成为 30。作者并没有在文章和 BioSample 网页中说明这一点，实则测序数据质控中质量和长度分布的部分可以省略。事实上，经过修建和质控的 SRA 数据有利于我们在有限的算力条件，使用较少的`fastq`组装出质量较高的完整基因组。我们仅使用覆盖度为 45.6x 的 SRR12072898 样本，即可得到全长 51,236,572 bp 的组装，与文章所给的 50,705,820 bp 误差为 0.60%。

在基因组组装过程中，得到了18747个 contigs，这些contigs高度碎片化，可能存在大量错误，特别是在预测总基因数方面。组装工具的参数设置会影响组装结果，推测是过于严格的参数导致更多的小contigs。此外，如果K-mer大小选择不当，也会影响组装的质量。一般来说，较长的K-mer可以更好地处理重复序列，但需要更高的计算资源，由于电脑内存不足选择了较短的K-mer。因此，在组装步骤中，可调整组装工具的参数，根据数据特点优化K-mer大小等设置；尝试使用不同的组装工具（Velvet、SOAPdenovo等），比较并优化结果。考虑到我们尝试复现论文中的实验方法，可根据论文中的标准，对组装结果进行长度过滤，只保留一定长度以上的contigs。使用去重复工具或算法，减少重复区域对组装的影响。

我们在搭建`RepeatMasker`软件和`RepeatModeler`软件的计算环境时，结合经验考虑到 RepeatModeler 依赖众多，选择单独构建虚拟环境，最终的软件包列表几乎没有冗余，在进行基因组组装中也沿用了这一规范。在进行其他步骤时，未充分考虑到不同步骤依赖的交叉范围，导致`pip`软件包过于繁杂，难以辨析不同步骤所需的环境配置，需要在以后的工作流程中加以改进。

重复序列和可移动元件预测结果显示，反转录元件93个，短散在核元件45个，长散在核元件25个，长末端重复序列(LTR) 23个。DNA转座子290个，简单重复序列2238个，低复杂度序列395个。从重复序列和可移动元件的分布来看，MC743 基因组中重复序列占比较低，仅为4.78%。其中，反转录元件和DNA转座子分别占基因组的0.29%和0.78%，简单重复序列和低复杂度序列分别占0.80%和0.15%。基因组中未分类的重复序列占比达2.68% / 1063个，提示存在大量尚未明确分类的重复序列。

研究者通过 MC743 与 MC848 CDSs 的比对，超过99%的基因分配到 10911 个唯一的正交群中，其中 6052个包含单拷贝正交基因，而剩下的 4859个包含来自同一菌株的至少2个基因。值得注意的是，1047个正交群在每个群中包含的基因数量上显示出菌株级别的差异。具体而言，在这1,047个正交群中，528个在MC848菌株中富集，497个在MC873菌株中富集，而剩下的 22个正交群包含菌株特异基因。MC848菌株显示了10个特异正交群，共包含20个基因，而其余12个MC873特异正交群共包含24个基因。我们的共线性分析显示，MC743 与 MC848 基因组在全区段具有显著的的正向互补匹配，这与两者相同物种的亲缘关系相符合。

我们起初采用构建BLAST数据库，将蛋白质序列与COG数据库比对的方法，实际运行中发现本地运行BLAST算法的效率过低，12 h后仅得到约20%的比对结果，且匹配度得分缺少过滤。因此，我们选择使用包含隐马尔可夫模型PFAM数据库进行 HMMER 比对，处理速率达到平均 3-4 seq/sec。这种方法的输出文件直观展示了蛋白质的全序列得分、最佳结构域得分和功能描述，便于数据过滤和处理，缺点在于 HMMER 功能注释缺少与 GO terms 的25种类型相似的明确分类，因此无法进行后续的 KEGG 代谢通路预测和蛋白质功能类型比对。

COG数据库包含了大量的蛋白质序列信息，总文件大小接近1GB。由于BLAST需要将整个数据库加载到内存中进行比对，如果内存不足就会导致运行效率非常低下。而且，BLAST算法本身是一种复杂的序列比对算法，对计算资源有较高的要求，其时间复杂度与序列长度和数据库大小成正比。当数据库非常大时，比对过程会变得十分耗时。针对这些问题，DIAMOND 是一种基于BLAST的序列比对算法，使用了一些优化技术，如基于种子的比对、双向最佳比对等，大幅提高比对效率。因此，在今后的学习中可以融合二者的优点，使用较快速的DIAMOND进行大致的比对，筛选出感兴趣的序列，再使用HMMER进行更细致的比对，兼顾效率和准确性。

### 参考文献

[1] Romeo O, Marchetta A, Giosa D, et al. Whole genome sequencing and comparative genome analysis of the halotolerant deep sea black yeast *Hortaea werneckii*[J]. Life, 2020, 10(10): 229. [Life | Free Full-Text | Whole Genome Sequencing and Comparative Genome Analysis of the Halotolerant Deep Sea Black Yeast Hortaea werneckii (mdpi.com)](https://www.mdpi.com/2075-1729/10/10/229)
[2] Shen Q, Chen Y, Jin D, et al. Comparative genome analysis of the oleaginous yeast *Trichosporon fermentans* reveals its potential applications in lipid accumulation[J]. Microbiological research, 2016, 192: 203-210. [Comparative genome analysis of the oleaginous yeast Trichosporon fermentans reveals its potential applications in lipid accumulation - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0944501316301793)
[3] Mo X, Zhou M, Li Y, et al. Safety assessment of a novel marine multi-stress-tolerant yeast *Meyerozyma guilliermondii* GXDK6 according to phenotype and whole genome-sequencing analysis[J]. Food Science and Human Wellness, 2024, 13(4): 2048-2059. [Safety assessment of a novel marine multi-stress-tolerant yeast Meyerozyma guilliermondii GXDK6 according to phenotype and whole genome-sequencing analysis - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2213453024000831)
[4] Baeza M, Zúñiga S, Peragallo V, et al. Identification of stress-related genes and a comparative analysis of the amino acid compositions of translated coding sequences based on draft genome sequences of Antarctic yeasts[J]. Frontiers in microbiology, 2021, 12: 623171. [Frontiers | Identification of Stress-Related Genes and a Comparative Analysis of the Amino Acid Compositions of Translated Coding Sequences Based on Draft Genome Sequences of Antarctic Yeasts (frontiersin.org)](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2021.623171/full)

### 附录

#### 计算环境配置信息

1. `repeatmasker`环境：用于第3步的重复序列和可移动元件预测。

   ```bash
   (repeatmasker) :$ pip list
   Package     Version
   ----------- ---------
   certifi     2022.12.7
   h5py        3.7.0
   mkl-fft     1.0.14
   mkl-random  1.0.4
   mkl-service 2.3.0
   numpy       1.17.0
   pip         22.3.1
   setuptools  65.6.3
   six         1.16.0
   wheel       0.38.4
   
   (repeatmasker) :$ conda list
   # packages in environment at /home/zfwang/anaconda3/envs/repeatmasker:
   #
   # Name                    Version                   Build  Channel
   _libgcc_mutex             0.1                        main    defaults
   _openmp_mutex             5.1                       1_gnu    defaults
   blas                      1.0                         mkl    defaults
   ca-certificates           2024.3.11            h06a4308_0    defaults
   certifi                   2022.12.7        py37h06a4308_0    defaults
   gdbm                      1.18                 hd4cb3f1_4    defaults
   h5py                      3.7.0            py37h737f45e_0    defaults
   hdf5                      1.10.6               hb1b8bf9_0    defaults
   intel-openmp              2023.1.0         hdb19cb5_46306    defaults
   ld_impl_linux-64          2.38                 h1181459_1    defaults
   libffi                    3.4.4                h6a678d5_0    defaults
   libgcc-ng                 11.2.0               h1234567_1    defaults
   libgfortran-ng            7.5.0               ha8ba4b0_17    defaults
   libgfortran4              7.5.0               ha8ba4b0_17    defaults
   libgomp                   11.2.0               h1234567_1    defaults
   libnsl                    2.0.0                h5eee18b_0    defaults
   libstdcxx-ng              11.2.0               h1234567_1    defaults
   mkl                       2019.4                      243    defaults
   mkl-service               2.3.0            py37he8ac12f_0    defaults
   mkl_fft                   1.0.14           py37hd81dba3_0    defaults
   mkl_random                1.0.4            py37hd81dba3_0    defaults
   ncurses                   6.4                  h6a678d5_0    defaults
   numpy                     1.17.0           py37h7e9f1db_0    defaults
   numpy-base                1.17.0           py37hde5b4d6_0    defaults
   openssl                   1.1.1w               h7f8727e_0    defaults
   perl                      5.26.2               h14c3975_0    defaults
   perl-app-cpanminus        1.7044                  pl526_1    bioconda
   perl-business-isbn        3.004                   pl526_0    bioconda
   perl-business-isbn-data   20140910.003            pl526_0    bioconda
   perl-carp                 1.38                    pl526_3    bioconda
   perl-common-sense         3.74                    pl526_2    bioconda
   perl-constant             1.33                    pl526_1    bioconda
   perl-data-dumper          2.173                   pl526_0    bioconda
   perl-devel-size           0.83            pl526h14c3975_0    bioconda
   perl-encode               2.88                    pl526_1    bioconda
   perl-exporter             5.72                    pl526_1    bioconda
   perl-extutils-makemaker   7.36                    pl526_1    bioconda
   perl-file-path            2.16                    pl526_0    bioconda
   perl-file-temp            0.2304                  pl526_2    bioconda
   perl-file-which           1.23                    pl526_0    bioconda
   perl-json                 4.02                    pl526_0    bioconda
   perl-json-xs              2.34                    pl526_1    bioconda
   perl-lwp                  6.35            pl526h14c3975_0    eumetsat
   perl-mime-base64          3.15                    pl526_1    bioconda
   perl-parent               0.236                   pl526_1    bioconda
   perl-storable             3.15            pl526h14c3975_0    bioconda
   perl-test-simple          1.302164                pl526_0    bioconda
   perl-uri                  1.76                    pl526_0    bioconda
   perl-xsloader             0.24                    pl526_0    bioconda
   pip                       22.3.1           py37h06a4308_0    defaults
   python                    3.7.16               h7a1cb2a_0    defaults
   readline                  8.2                  h5eee18b_0    defaults
   setuptools                65.6.3           py37h06a4308_0    defaults
   six                       1.16.0             pyhd3eb1b0_1    defaults
   sqlite                    3.41.2               h5eee18b_0    defaults
   tk                        8.6.12               h1ccaba5_0    defaults
   wheel                     0.38.4           py37h06a4308_0    defaults
   xz                        5.4.6                h5eee18b_0    defaults
   zlib                      1.2.13               h5eee18b_0    defaults
   ```

2. `bio2503`环境：用于第1步、第4-7步的质控和基因组分析。

   ```bash
   (bio2503) :$ pip list
   Package                           Version
   --------------------------------- ------------
   anaconda-anon-usage               0.4.4
   anaconda-catalogs                 0.2.0
   anaconda-client                   1.12.3
   anaconda-cloud-auth               0.1.4
   anaconda-navigator                2.5.2
   anaconda-project                  0.11.1
   biopython                         1.83
   cachetools                        4.2.2
   certifi                           2024.2.2
   conda                             23.1.0
   conda-build                       24.1.2
   conda-content-trust               0.2.0
   conda_index                       0.4.0
   conda-libmamba-solver             22.8.1
   conda-pack                        0.6.0
   conda-package-handling            2.2.0
   conda_package_streaming           0.9.0
   conda-repo-cli                    1.0.75
   conda-token                       0.4.0
   conda-verify                      3.4.2
   dask                              2023.11.0
   datashader                        0.16.0
   debugpy                           1.6.7
   distlib                           0.3.8
   distributed                       2023.11.0
   distro                            1.8.0
   docstring-to-markdown             0.11
   docutils                          0.18.1
   entrypoints                       0.4
   et-xmlfile                        1.1.0
   executing                         0.8.3
   fastjsonschema                    2.16.2
   filelock                          3.13.1
   flake8                            6.0.0
   Flask                             2.2.5
   fonttools                         4.25.0
   frozenlist                        1.4.0
   fsspec                            2023.10.0
   future                            0.18.3
   gensim                            4.3.0
   gitdb                             4.0.7
   GitPython                         3.1.37
   gmpy2                             2.1.2
   greenlet                          3.0.1
   h5py                              3.9.0
   HeapDict                          1.0.1
   holoviews                         1.18.3
   hvplot                            0.9.2
   hyperlink                         21.0.0
   ipykernel                         6.28.0
   ipython                           8.20.0
   ipython-genutils                  0.2.0
   ipywidgets                        7.6.5
   itemadapter                       0.3.0
   itemloaders                       1.1.0
   itsdangerous                      2.0.1
   jupyter                           1.0.0
   jupyter_client                    8.6.0
   jupyter-console                   6.6.3
   jupyter_core                      5.5.0
   jupyter-events                    0.8.0
   jupyter-lsp                       2.2.0
   jupyter_server                    2.10.0
   jupyter_server_terminals          0.4.4
   jupyterlab                        4.0.11
   jupyterlab-pygments               0.1.2
   jupyterlab_server                 2.25.1
   jupyterlab-widgets                3.0.9
   keyring                           23.13.1
   kiwisolver                        1.4.4
   lazy_loader                       0.3
   lazy-object-proxy                 1.6.0
   lckr_jupyterlab_variableinspector 3.1.0
   libarchive-c                      2.9
   libmambapy                        1.1.0
   Markdown                          3.4.1
   markdown-it-py                    2.2.0
   MarkupSafe                        2.1.3
   matplotlib                        3.8.4
   matplotlib-inline                 0.1.6
   mkl-fft                           1.3.8
   mkl-random                        1.2.4
   mkl-service                       2.4.0
   mypy                              1.8.0
   mypy-extensions                   1.0.0
   navigator-updater                 0.4.0
   notebook                          7.0.8
   notebook_shim                     0.2.3
   numba                             0.59.1
   numexpr                           2.8.7
   numpy                             1.26.4
   numpydoc                          1.5.0
   openpyxl                          3.0.10
   overrides                         7.4.0
   packaging                         23.1
   pandas                            2.1.4
   pip                               23.3.1
   pipenv                            2023.12.1
   pkce                              1.0.3
   pkginfo                           1.9.6
   pycurl                            7.45.2
   python-dateutil                   2.8.2
   python-dotenv                     0.21.0
   python-json-logger                2.0.7
   python-lsp-black                  2.0.0
   python-lsp-jsonrpc                1.1.2
   python-lsp-server                 1.7.2
   python-slugify                    5.0.2
   python-snappy                     0.6.1
   pytoolconfig                      1.2.6
   queuelib                          1.6.2
   referencing                       0.30.2
   regex                             2023.10.3
   requests                          2.31.0
   requests-file                     1.5.1
   requests-toolbelt                 1.0.0
   rfc3339-validator                 0.1.4
   rfc3986-validator                 0.1.1
   Rtree                             1.0.1
   ruamel.yaml                       0.17.21
   ruamel-yaml-conda                 0.17.21
   service-identity                  18.1.0
   setuptools                        68.2.2
   six                               1.16.0
   snowballstemmer                   2.2.0
   sortedcontainers                  2.4.0
   soupsieve                         2.5
   Sphinx                            5.0.2
   sphinxcontrib-applehelp           1.0.2
   sphinxcontrib-devhelp             1.0.2
   sphinxcontrib-htmlhelp            2.0.0
   sphinxcontrib-jsmath              1.0.1
   sphinxcontrib-qthelp              1.0.3
   sphinxcontrib-serializinghtml     1.1.5
   spyder                            5.4.3
   spyder-kernels                    2.4.4
   sympy                             1.12
   tables                            3.9.2
   tabulate                          0.9.0
   toolz                             0.12.0
   tornado                           6.3.3
   Twisted                           23.10.0
   typing_extensions                 4.9.0
   tzlocal                           2.1
   uc-micro-py                       1.0.1
   ujson                             5.4.0
   Unidecode                         1.2.0
   urllib3                           2.0.7
   validators                        0.18.2
   wheel                             0.41.2
   widgetsnbextension                3.5.2
   wrapt                             1.14.1
   wurlitzer                         3.0.2
   yapf                              0.31.0
   yarl                              1.9.3
   zict                              3.0.0
   zipp                              3.17.0
   
   (bio2503) :$ conda list
   # packages in environment at /home/zfwang/anaconda3/envs/bio2503:
   #
   # Name                    Version                   Build  Channel
   _libgcc_mutex             0.1                        main    defaults
   _openmp_mutex             5.1                       1_gnu    defaults
   bedtools                  2.26.0                        0    bioconda
   blast                     2.5.0                hc0b0e79_3    bioconda
   boost                     1.57.0                        4    https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   bzip2                     1.0.6                         3    https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   ca-certificates           2024.3.11            h06a4308_0    defaults
   curl                      7.26.0                        1    https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   hmmer                     3.1b2                         3    bioconda
   icu                       54.1                          0    https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   libdeflate                1.2                  h516909a_1    bioconda
   libgcc                    7.2.0                h69d50b8_2    defaults
   libgcc-ng                 11.2.0               h1234567_1    defaults
   libgomp                   11.2.0               h1234567_1    defaults
   libstdcxx-ng              11.2.0               h1234567_1    defaults
   ncurses                   5.9                          10    https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   openssl                   1.1.1w               h7f8727e_0    defaults
   samtools                  1.8                           4    bioconda
   xz                        5.2.3                         0    https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   zlib                      1.2.11                        0    https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
   ```

3. `quast_env`环境：用于第2步基因组组装。

   ```bash
   (quast_env) :$ pip list
   Package           Version
   ----------------- --------
   certifi           2024.2.2
   cycler            0.11.0
   fonttools         4.38.0
   kiwisolver        1.4.2
   matplotlib        3.5.3
   numpy             1.21.6
   packaging         23.2
   Pillow            9.2.0
   pip               22.3.1
   pyparsing         3.1.2
   python-dateutil   2.9.0
   setuptools        65.6.3
   six               1.16.0
   tornado           6.2
   typing_extensions 4.7.1
   unicodedata2      14.0.0
   wheel             0.38.4
   
   (quast_env) :$ conda list
   # packages in environment at /home/zfwang/anaconda3/envs/bio2503:
   #
   # Name                    Version                   Build  Channel
   _libgcc_mutex             0.1                 conda_forge    conda-forge
   _openmp_mutex             4.5                       2_gnu    conda-forge
   brotli                    1.1.0                hd590300_1    conda-forge
   brotli-bin                1.1.0                hd590300_1    conda-forge
   bzip2                     1.0.8                hd590300_5    conda-forge
   ca-certificates           2024.3.11            h06a4308_0    defaults
   certifi                   2024.2.2           pyhd8ed1ab_0    conda-forge
   cycler                    0.11.0             pyhd8ed1ab_0    conda-forge
   dbus                      1.13.6               h5008d03_3    conda-forge
   expat                     2.6.2                h59595ed_0    conda-forge
   fontconfig                2.14.2               h14ed4e7_0    conda-forge
   fonttools                 4.38.0           py37h540881e_0    conda-forge
   freetype                  2.12.1               h267a509_2    conda-forge
   glib                      2.80.2               hf974151_0    conda-forge
   glib-tools                2.80.2               hb6ce0ca_0    conda-forge
   gst-plugins-base          1.14.0               hbbd80ab_1    defaults
   gstreamer                 1.14.1               h5eee18b_1    defaults
   icu                       58.2              hf484d3e_1000    conda-forge
   jpeg                      9e                   h0b41bf4_3    conda-forge
   kiwisolver                1.4.2            py37h7cecad7_1    conda-forge
   lcms2                     2.14                 h6ed2654_0    conda-forge
   ld_impl_linux-64          2.38                 h1181459_1    defaults
   lerc                      3.0                  h9c3ff4c_0    conda-forge
   libblas                   3.9.0           22_linux64_openblas    conda-forge
   libbrotlicommon           1.1.0                hd590300_1    conda-forge
   libbrotlidec              1.1.0                hd590300_1    conda-forge
   libbrotlienc              1.1.0                hd590300_1    conda-forge
   libcblas                  3.9.0           22_linux64_openb
   ```

#### 项目核心脚本

以下所有脚本和项目报告 Markdown 已发布到 [github个人主页](https://github.com/BleauRime0/Deep-Sea-Black-Yeast)。

1. `fastq_analysis.sh`：用于解析 SRA 文件的大小和短读段长度。

   ```bash
   #!/bin/bash
   
   DIRECTORY="~/biolab/quality_control/"
   
   analyze_fastq() {
       for file in "$DIRECTORY"/*.fastq; do
           if [ -f "$file" ]; then
               echo "Analyzing file: $(basename "$file")"
               analyze_file "$file"
           fi
       done
   }
   
   analyze_file() {
       local file="$1"
       local total_reads=0
       local min_length=999999
       local max_length=0
   
       while read -r line; do
           if [[ "$((++total_reads % 4))" == "2" ]]; then
               local length=${#line}
               if [[ "$length" -lt "$min_length" ]]; then
                   min_length="$length"
               fi
               if [[ "$length" -gt "$max_length" ]]; then
                   max_length="$length"
               fi
           fi
       done < "$file"
   
       echo "Total Reads: $total_reads"
       echo "Read Length: $min_length - $max_length"
   }
   
   analyze_fastq
   ```

2. `extract_scaffold.py`：在蛋白质编码序列预测中，对`GTF`和`fna`文件的注释部分进行处理，保留`scaffold`编号作为识别特征，适用于`gffread`工具。

   ```python
   def description_process(input_file, output_file):
       with open(input_file, 'r') as input_file, open(output_file, 'w') as output_file:
           for line in input_file:
               if line.startswith('>'):
                   parts = line.strip().split(',')
                   output_file.write(f">{parts[0][-13:]}\n")
               else:
                   output_file.write(line)
   
   def genemark_process(input_file, output_file):
       with open(input_file, 'r') as input_file, open(output_file, 'w') as output_file:
           for line in input_file:
               fields = line.strip().split('\t')
               first_field = fields[0]
               comma_index = first_field.find(',')
               new_first_field = first_field[:comma_index]
               if len(new_first_field) > 13:
                   new_first_field = new_first_field[-13:]
               tab_separated_fields = "\t".join(fields[1:])
               output_line = f"{new_first_field}\t{tab_separated_fields}"
               output_file.write(output_line + '\n')
   
   description_process('MC743.fna', 'MC743.fsa')
   genemark_process('MC743.genemark.gtf', 'MC743.gtf')
   ```

3. `hmmer_results.py`：在蛋白质功能预测中，对`HMMER`比对结果文件`_COG_hmmer.txt`进行过滤和特征提取，输出蛋白质功能表格。

   ```python
   import csv
   
   def filter_hmmer_output(input_file, output_file, evalue_threshold=1e-10, score_threshold=0):
       with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
           tsv_writer = csv.writer(outfile, delimiter='\t')
           tsv_writer.writerow(['Gene Name', 'Accession', 'Protein Sequence ID', 'E value', 'Score', 'Bias', 'Protein Function'])
   
           for line in infile:
               if line.startswith('#') or not line.strip():
                   continue
   
               parts = line.split()
               target_name = parts[0]
               accession = parts[1]
               query_name = parts[2]
               evalue = float(parts[4])
               score = float(parts[5])
               bias = float(parts[6])
               description = ' '.join(parts[18:])
   
               if evalue <= evalue_threshold and score-bias >= score_threshold:
                   tsv_writer.writerow([target_name, accession, query_name, evalue, score, bias, description])
   
   filter_hmmer_output('MC743_COGs_hmmer.txt', 'MC743_protein_function.tsv')
   ```
