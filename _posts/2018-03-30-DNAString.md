---
title: "bioconductor and markdown"
date: 2018-03-30 14:56:09
categories:
- R
- bioconductor
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### bioconductor

#### **下载**  
```{r cars,eval=FALSE}
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicFeatures", "GenomicRanges"))
``` 
      
##### XString是一个“虚拟类”，不能被“实例化”（不能创建XString对象）。四个子类(对象操作类似)为：
 - BString：用于存储一般字符串。
 - DNAString：用于存储DNA（核苷酸）序列。
 - RNAString：用于存储RNA序列。
 - AAString：用于存储蛋白质（氨基酸）序列。
      
#### **创建Biostrings的对象:**
```{r results="hide",message=FALSE}   
library(Biostrings)
```
```{r}
a=BString("I am a string!")
a
length(a) 
a[1:4]          #子集的取出Subsetting:
subseq(a,1,4)   #效率优于a[1:4]
```  
     
#### **可以将序列倒转:**
```{r}
rev(a)
```  
     
#### **比较判断以及转存到字符串(实)中:**
```{r results='markup'}
a=="I am"
a[1:4]=="I am"
toString(a)
class(a)
class(toString(a))
```

#### **DNAString/RNAString**
##### 他们只能采取“有效”字符来表示核苷酸：
```{r}
IUPAC_CODE_MAP   #显示核苷酸与表示核苷酸的有效字符的一一对应关系
DNA_ALPHABET
DNA_BASES
RNA_ALPHABET
RNA_BASES
```  
  
a=DNAString(“I am a string”)会报错，因为I等不属于默认有效字符    
b=RNAString(“ATTGCC”)也会报错，因为RNA_BASES中默认不含T    
     
```{r}
a=DNAString("ATTGCC")
a
b=RNAString("AUUGCC")
b
complement(a)   #互补链
reverse(a)      #反向链
reverseComplement(a)      #反向互补链
translate(a)    #翻译序列，以氨基酸缩写显示
```

#### **碱基频率计数**
```{r}
alphabetFrequency(a)
alphabetFrequency(a, baseOnly=TRUE)    #屏蔽MN等（M.N...是什么？）
letterFrequency(a,"C")
letterFrequency(a,"CGT")
```

#### **单字符串匹配和对齐**
##### 功能分为以下四组：   
 - 查找给定模式的出现次数：matchPattern，countPattern，vmatchPattern，vcountPattern   
 - 根据参考匹配模式字典：matchPDict，countPDict   
 - 与位置匹配/计数权重矩阵（PWM）：matchPWM，countPWM，PWMscoreStartingAt   
 - 全局/局部比对：pairwiseAlignment，stringDist   

#### **matchPattern**
##### 在序列中找寻一个给定的模式，可允许失配和插入/缺失（插入缺失）：
```{r}
a=DNAString("ACGTACGTACGC")
matchPattern("CGT", a)
matchPattern("CGT", a, max.mismatch=1)    #允许最大错配为1
m=matchPattern("CGT", a, max.mismatch=1)
start(m)    #每个匹配到的mers的起始点
end(m)      #每个匹配到的mers的末尾
length(m)     #有多少个匹配到该模式的mers
countPattern("CGT", a, max.mismatch=1)    #有多少个匹配到该模式的mers
```
##### 这些功能可以有效地计算庞大基因组中的n-mers存在，例如：   
 - GC含量：“C”出现数+“G”出现数（当然利用之前所说的letterFrequency频数功能更有效）   
 - CpG含量：发生“CG”的数量    

#### **matchPDict**
##### 查找出现的一组模式。或者你可以写一个循环，但matchPDict更高效（R循环是出了名慢）。
```{r}
a=DNAString("ACGTACGTACGC")
dict0=PDict(c("CGT","ACG"))
mm=matchPDict(dict0, a)
mm[[1]]   #匹配CGT的情况
mm[[2]]   #匹配ACG的情况
```   
    
####  **PWM**
#####  PWM(位置权重矩阵), 用来表示DNA motif
 - motif是蛋白中较小的保守序列片断   
```{r}
a=DNAString("ACGTACGTACTC")
motif=matrix(c(0.97,0.01,0.01,0.01,0.1,0.5,0.39,0.01,0.01,0.05,0.5,0.44),nrow=4)  #建立得分矩阵
rownames(motif)=c("A","C","G","T")     #得分矩阵行名（A在第1、2、3个位子的分数...所以默认是3个碱基的序列片段）
motif
matchPWM(motif, a)
countPWM(motif, a)   #找到的motif个数
PWMscoreStartingAt(motif, a, 1:10)
```
#### **多个字符串（序列）的操作String views and set**
 - 可以在一个循环中实现多个字符串的操作，但效率非常低。   
 - 多个字符串来源于“mother”字符串，并且放入一个字符串“view”或“set”。   
 - XStringViews：包含多个“view”（相同字符串的start/end位置）。   
 - DNAStringSet/RNAStringSet：类似但创建实际DNA/RNAString实例。   
 - StringSet允许比StringViews更多的操作。      

#### **XStringViews的基本操作**
```{r}
a=DNAString("ACGTACGTACTC")
a2=Views(a, start=c(1,5,8), end=c(3,8,12))
a2
subject(a2)    #a2的来源主体a
length(a2)    #a2有几个
start(a2)     #a2分别的开始位置
end(a2)       #a2分别的结束位置
alphabetFrequency(a2, baseOnly=TRUE)    #分别统计a2中每个字符串中的碱基频率
a2==DNAString("ACGT")    #判断是否等于（逻辑）
toString(a2)              #提出已设定位置的每个序列a2
```
#### **DNAStringSet的基本操作**
```{r}
a=DNAString("ACGTACGTACTC")
a2=DNAStringSet(a, start=c(1,5,9), end=c(4,8,12))
a2
a2[[3]]
alphabetFrequency(a2, baseOnly=TRUE)
```

#### **一些操作只在StringSet中被允许，而不存在于Views，如set操作**
```{r error=TRUE}
a1=DNAStringSet(a, start=c(1,5,9), end=c(4,8,12))
a1
unique(a1)      #去掉了提取出来的重复的片段
a2=Views(a, start=c(1,5,9), end=c(4,8,12))    
unique(a2)     #该句会报错，Views不允许unique操作
a1=Views(a, start=c(1,9), end=c(4,12))
a1
a2=Views(a, start=c(1), end=c(4))
a2
setdiff(a1,a2)   #报错
union(a1, a2)    #报错
a1=DNAStringSet(a, start=c(1,9), end=c(4,12))
a1
a2=DNAStringSet(a, start=c(1), end=c(4))
a2
setdiff(a1,a2)    #差集
union(a1,a2)    #并集
```
#### **多序列（字符串）的匹配** 
 - 使用vmatchPattern和vmatchPDict   
 - 对于PWM没有对应的功能   
```{r error=TRUE}
a=DNAString("ACGTACGTACTC")
a2=DNAStringSet(a, start=c(1,5,9), end=c(4,8,12))
a2
vv=vmatchPattern("CG", a2)
vv      #对于a2的每一行都匹配搜索CG，如果某一行没有CG，则不显示数字；显示232表示a2中这一行的序列CG在第2、3个位置
vv[[1]]    
# 这些同样不能在Views中操作   
a2=Views(a, start=c(1,5,9), end=c(4,8,12))
a2
vv=vmatchPattern("CG", a2)
```