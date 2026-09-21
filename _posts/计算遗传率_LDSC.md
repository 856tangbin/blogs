# LDSC 学习与实践笔记

**作者**：吴玥冉
**日期**：2026年3月
**工具**：R语言以及cmd命令符

------

## 前言

LDSC（LD Score Regression）是进行遗传率估算和遗传相关分析的重要工具。本文档记录了我学习 LDSC 的过程，包括原理理解、环境配置、数据准备、常见问题以及最终结果解读。其中包含两个性状的实际分析流程：微生物组（c_Bacilli）和代谢组（1-methyl-histidine）。

------

##  一、LDSC 原理简介

LDSC 利用连锁不平衡（LD）结构，通过回归 SNP 的卡方统计量（χ²）对 LD Score 的关系，来估计性状的遗传率（h²）和判断群体分层。其核心公式为：

E[χj2]=1+N⋅h2M⋅ℓj+N⋅τ*E*[*χ**j*2]=1+*N*⋅*M**h*2⋅ℓ*j*+*N*⋅*τ*

- χj2*χ**j*2：SNP j 的卡方统计量
- N*N*：样本量
- h2*h*2：遗传率
- M*M*：SNP 总数
- ℓjℓ*j*：SNP j 的 LD Score（与周围 SNP 的 r² 之和）
- τ*τ*：控制其他偏误（如群体分层）

**关键点**：

- 斜率 → 遗传率
- 截距 → 群体分层（越接近 1 越好）

------

## 二、安装与配置

### 2.1 创建 Python 2.7 环境

LDSC 官方代码基于 Python 2.7，使用 conda 创建独立环境。

bash

```
conda create -n ldsc python=2.7
conda activate ldsc
```



### 2.2 安装 LDSC

从 GitHub 下载代码：

bash

```
git clone https://github.com/bulik/ldsc.git
cd ldsc
```



安装依赖包：

bash

```
pip install -r requirements.txt
```



估计欧洲千人基因组（1000 Genomes Europeans）的HapMap3 LD Scores: 







### 2.3 下载参考文件

LDSC 需要参考 LD 分数文件和 HapMap3 SNP 列表。

- **东亚人群 LD 分数**：`eas_w_ld_chr.tar.bz2`
- **HapMap3 SNP 列表**：`w_hm3.snplist.bz2`



解压后，将文件放置于 `ldsc/ref/eas_w_ld_chr/` 和 `ldsc/ref/w_hm3.snplist`。

------

##  三、数据准备（微生物组）

### 3.1 坐标转 rsID

使用 R 包 `SNPlocs.Hsapiens.dbSNP144.GRCh37` 进行转换（基于 GRCh37，dbSNP 144 版本）。

**R 代码**：

r

```
library(data.table)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

setwd("D:/桌面/ldsc")

gwas <- fread("c_Bacilli_clean.csv")
gwas[, c("chr", "pos") := tstrsplit(SNP, ":", fixed = TRUE)]
gwas[, chr := gsub("chr", "", chr)]
gwas[, pos := as.integer(pos)]

chrs <- intersect(unique(gwas$chr), c(1:22, "X", "Y"))
tolerance <- 100
matched_list <- list()

for (chr_i in chrs) {
  cat("\nProcessing chr", chr_i)
  gwas_chr <- gwas[chr == chr_i, .(pos)]
  if (nrow(gwas_chr) == 0) next
  chr_snps <- snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, as.character(chr_i))
  dt_snps <- as.data.table(chr_snps)[, .(pos, rsid = RefSNP_id)]
  
  matched <- dt_snps[gwas_chr, .(pos = i.pos, rsid = x.rsid), 
                     on = .(pos >= i.pos - tolerance, pos <= i.pos + tolerance), 
                     allow.cartesian = FALSE, mult = "first"]
  matched <- matched[!is.na(rsid)]
  if (nrow(matched) > 0) {
    matched[, chr := chr_i]
    matched_list[[as.character(chr_i)]] <- matched
  }
}
matched_all <- rbindlist(matched_list)
gwas_final <- merge(gwas, matched_all, by = c("chr", "pos"))
gwas_final[, SNP := rsid]
gwas_final[, c("chr", "pos", "rsid") := NULL]
fwrite(gwas_final, "c_Bacilli_with_rsid_final.csv")
```



运行后，成功匹配 **9,734,986** 个 SNP。

### 3.2 与参考面板合并（验证过程）

r

```
ref_dir <- "D:/桌面/ldsc/ref/eas_w_ld_chr/"
ld_files <- list.files(ref_dir, pattern = "\\.l2\\.ldscore\\.gz$", full.names = TRUE)
ref_snps <- unique(rbindlist(lapply(ld_files, fread, select = "SNP")))
cat("参考面板 rsID 数：", nrow(ref_snps), "\n")

gwas_rs <- fread("c_Bacilli_with_rsid_final.csv")
gwas_ready <- merge(gwas_rs, ref_snps, by = "SNP", all = FALSE)
cat("合并后 SNP 数：", nrow(gwas_ready), "\n")
fwrite(gwas_ready, "c_Bacilli_ready_for_munge.csv")
```



结果显示，通过简单 rsID 匹配可获得 **100,292** 个重叠 SNP。

### 3.3 生成 `.sumstats.gz`（Python 脚本）

由于 `munge_sumstats.py` 一直报错，改用 Python 脚本手动生成。（由deepseek生成）

python

```
import gzip

input_file = 'c_Bacilli_ready_for_munge.csv'
output_file = 'c_Bacilli_final.sumstats.gz'

with open(input_file, 'r') as f_in, gzip.open(output_file, 'wb') as f_out:
    f_in.readline()  # skip header
    f_out.write(b'SNP\tA1\tA2\tZ\tP\tN\n')
    for line in f_in:
        parts = line.strip().split(',')
        snp = parts[0]
        a1 = parts[1]
        a2 = parts[2]
        beta = float(parts[4])
        se = float(parts[5])
        p = float(parts[6])
        n = float(parts[7])
        z = beta / se
        out_line = f"{snp}\t{a1}\t{a2}\t{z}\t{p}\t{n}\n"
        f_out.write(out_line.encode('utf-8'))
print("Done")
```



------

##  四、数据准备

### 4.1 坐标转 rsID（完整代码）

代谢组与微生物组方法完全一致，仅输入输出文件名不同。以下为完整代码：

r

```
library(data.table)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

setwd("D:/桌面/ldsc")

metabo <- fread("1-methyl-histidine_clean.csv")
metabo[, c("chr", "pos") := tstrsplit(SNP, ":", fixed = TRUE)]
metabo[, chr := gsub("chr", "", chr)]
metabo[, pos := as.integer(pos)]

chrs <- intersect(unique(metabo$chr), c(1:22, "X", "Y"))
tolerance <- 100
matched_list <- list()

for (chr_i in chrs) {
  cat("\nProcessing chr", chr_i)
  gwas_chr <- metabo[chr == chr_i, .(pos)]
  if (nrow(gwas_chr) == 0) next
  chr_snps <- snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37, as.character(chr_i))
  dt_snps <- as.data.table(chr_snps)[, .(pos, rsid = RefSNP_id)]
  
  matched <- dt_snps[gwas_chr, .(pos = i.pos, rsid = x.rsid), 
                     on = .(pos >= i.pos - tolerance, pos <= i.pos + tolerance), 
                     allow.cartesian = FALSE, mult = "first"]
  matched <- matched[!is.na(rsid)]
  if (nrow(matched) > 0) {
    matched[, chr := chr_i]
    matched_list[[as.character(chr_i)]] <- matched
  }
}
matched_all <- rbindlist(matched_list)
metabo_final <- merge(metabo, matched_all, by = c("chr", "pos"))
metabo_final[, SNP := rsid]
metabo_final[, c("chr", "pos", "rsid") := NULL]
fwrite(metabo_final, "1methyl_with_rsid_final.csv")
```



### 4.2 与参考面板合并（验证）

r

```
ref_dir <- "D:/桌面/ldsc/ref/eas_w_ld_chr/"
ld_files <- list.files(ref_dir, pattern = "\\.l2\\.ldscore\\.gz$", full.names = TRUE)
ref_snps <- unique(rbindlist(lapply(ld_files, fread, select = "SNP")))
cat("参考面板 rsID 数：", nrow(ref_snps), "\n")

metabo_rs <- fread("1methyl_with_rsid_final.csv")
metabo_ready <- merge(metabo_rs, ref_snps, by = "SNP", all = FALSE)
cat("合并后 SNP 数：", nrow(metabo_ready), "\n")
fwrite(metabo_ready, "1methylhistidine_ready_for_munge.csv")
```



同样得到 **100,292** 个重叠 SNP。

### 4.3 生成 `.sumstats.gz`（Python 脚本）

python

```
import gzip

input_file = '1methylhistidine_ready_for_munge.csv'
output_file = '1methylhistidine_final.sumstats.gz'

with open(input_file, 'r') as f_in, gzip.open(output_file, 'wb') as f_out:
    f_in.readline()  # skip header
    f_out.write(b'SNP\tA1\tA2\tZ\tP\tN\n')
    for line in f_in:
        parts = line.strip().split(',')
        snp = parts[0]
        a1 = parts[1]
        a2 = parts[2]
        beta = float(parts[4])
        se = float(parts[5])
        p = float(parts[6])
        n = float(parts[7])
        z = beta / se
        out_line = f"{snp}\t{a1}\t{a2}\t{z}\t{p}\t{n}\n"
        f_out.write(out_line.encode('utf-8'))
print("Done")
```



------

##  五、运行 LDSC

### 5.1 通用模板命令

LDSC 运行的基本命令格式如下：

bash

```
python ldsc.py --h2 <sumstats_file> --ref-ld-chr <ref_dir> --w-ld-chr <ref_dir> --out <output_prefix>
```



### 5.2 微生物组

**生成 sumstats**（执行 Python 脚本）：

bash

```
conda activate ldsc
cd D:\桌面\ldsc
python make_sumstats_microbe.py   # 对应上述微生物组 Python 脚本，文件名为 make_sumstats_microbe.py
```



**运行 LDSC**：

bash

```
python ldsc.py --h2 c_Bacilli_final.sumstats.gz --ref-ld-chr ref/eas_w_ld_chr/ --w-ld-chr ref/eas_w_ld_chr/ --out c_Bacilli_final_h2
```



### 5.3 代谢组

**生成 sumstats**（执行 Python 脚本）：

bash

```
conda activate ldsc
cd D:\桌面\ldsc
python make_sumstats_metab.py   # 对应上述代谢组 Python 脚本，文件名为 make_sumstats_metab.py
```



**运行 LDSC**：

bash

```
python ldsc.py --h2 1methylhistidine_final.sumstats.gz --ref-ld-chr ref/eas_w_ld_chr/ --w-ld-chr ref/eas_w_ld_chr/ --out 1methylhistidine_final_h2
```



------

## 六、结果解读

### 6.1 微生物组结果

text

```
After merging with reference panel LD, 88011 SNPs remain.
Total Observed scale h2: -0.481 (0.3956)
Lambda GC: 0.9976
Mean Chi^2: 1.0021
Intercept: 1.0172 (0.0139)
```



### 6.2 代谢组结果

text

```
After merging with reference panel LD, 88011 SNPs remain.
Total Observed scale h2: -0.6846 (0.3581)
Lambda GC: 1.0074
Mean Chi^2: 1.0018
Intercept: 1.0251 (0.0136)
```



### 6.3 对比分析

| 指标               | 微生物组        | 代谢组           |
| :----------------- | :-------------- | :--------------- |
| 合并后 SNP 数      | 88,011          | 88,011           |
| 遗传率 h² (标准误) | -0.481 (0.3956) | -0.6846 (0.3581) |
| 平均卡方           | 1.0021          | 1.0018           |
| 截距               | 1.0172 (0.0139) | 1.0251 (0.0136)  |
| λGC                | 0.9976          | 1.0074           |

分析：

1.两个性状均未表现出可检测的遗传基础。平均卡方接近 1，表明 GWAS 统计量中缺乏多基因信号。

2.遗传率估计为负值，说明结果不可靠，但负值本身提示真实遗传率可能接近 0。

3.截距略高于 1，提示可能存在轻微群体分层，但不影响主要结论。

4.数据质量（λGC、截距）尚可，但样本量可能不足以捕捉微弱信号，或性状本身确实不受遗传因素影响。

------

##  七、常见问题与解决方案

| 问题                              | 解决方案                                                     |
| --------------------------------- | :----------------------------------------------------------- |
| R 包安装失败（依赖问题）          | 手动下载 `.tar.gz` 文件，用 `install.packages` 本地安装；或从 CRAN 安装依赖 |
| `munge_sumstats.py` 找不到 SNP 列 | 改用 Python 脚本手动生成 sumstats.gz（按列索引读取）         |
| 坐标转 rsID 匹配率低              | 使用范围匹配（±100bp），避免因位置微小偏移漏掉匹配           |

------

##  八、实践总结

1.生信要多实操：光看教程很容易忽略细节，只有亲手跑一遍才能遇到真实问题。

2.不要死磕 munge：如果 munge 反复报列名错误，经各种各样的排查无果，最后求助于ai用 Python 手动处理。

3.范围匹配很有用：位置误差允许 ±100bp 能显著提高匹配率，尤其是不同版本或软件间微小偏移。增加范围后匹配数增加了很多，从一万六千提升到了十万左右。

4.R 包安装需耐心：R包的安装有时候也需要特定的环境

5.验证合并数：通过 R 直接比较 rsID 文件与参考面板的 SNP 列，可以提前预估 LDSC 的可用 SNP 数量，并理解手动匹配与最终使用数量的差异。

