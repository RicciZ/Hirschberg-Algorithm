# CS466 Introduction to Bioinformatics Course Project - Analysis and Implementation of Hirschberg Algorithm

## Contents
  - [Abstract](#abstract)
  - [Global Alignment](#global-alignment)
  - [Needleman-Wunsch](#needleman-wunsch)
  - [Hirschberg](#hirschberg)
  - [Analysis](#analysis)
  - [Run](#run)


<a name="abstract"></a>
## Abstract
In this course project for CS466 Introduction to Bioinformatics, we implemented and analyzed the time complexity and space complexity of the linear space global alignment algorithm - Hirschberg algorithm, and made a comparison to the Needleman-Wunsch algorithm. 
<pre xml:lang="latex">\sqrt{2}</pre>
$'\sqrt{2}'$
```math
\sqrt{2}
```
<a name="global"></a>
## Global Alignment
Given two strings <img src="https://latex.codecogs.com/svg.image?\boldsymbol{v}&space;\in&space;\Sigma^m" title="\boldsymbol{v} \in \Sigma^m" /> and <img src="https://latex.codecogs.com/svg.image?\boldsymbol{w}&space;\in&space;\Sigma^n" title="\boldsymbol{w} \in \Sigma^n" />, where <img src="https://latex.codecogs.com/svg.image?\Sigma" title="\Sigma" /> is the set of characters, a scoring function <img src="https://latex.codecogs.com/svg.image?\delta:(\Sigma&space;\cup&space;\{-\})\times(\Sigma&space;\cup&space;\{-\})&space;\rightarrow&space;\mathbb{R}" title="\delta:(\Sigma \cup \{-\})\times(\Sigma \cup \{-\}) \rightarrow \mathbb{R}" />, find alignment with maximum score.

An alignment <img src="https://latex.codecogs.com/svg.image?\boldsymbol{A}=[a_{i,j}]" title="\boldsymbol{A}=[a_{i,j}]" /> is a <img src="https://latex.codecogs.com/svg.image?2&space;\times&space;k" title="2 \times k" /> matrix s.t. 
- <img src="https://latex.codecogs.com/svg.image?k&space;=&space;\{max(m,n),\cdots,m&plus;n\}" title="k = \{max(m,n),\cdots,m+n\}" />
- <img src="https://latex.codecogs.com/svg.image?a_{i,j}&space;\in&space;\Sigma&space;\cup&space;\{-\}" title="a_{i,j} \in \Sigma \cup \{-\}" />
- there is no <img src="https://latex.codecogs.com/svg.image?j&space;\in&space;[k]" title="j \in [k]" /> where <img src="https://latex.codecogs.com/svg.image?a_{1,j}=a_{2,j}=-" title="a_{1,j}=a_{2,j}=-" />

<a name="needleman-wunsch"></a>
## Needleman-Wunsch
<img src="https://latex.codecogs.com/svg.image?s[i,j]&space;=&space;\max\begin{cases}&space;&space;0&space;&\text{if&space;}&space;i=0&space;\text{&space;and&space;}&space;j=0&space;\\&space;&space;s[i-1,j]&space;&plus;&space;\delta(v_i,-)&space;&\text{if&space;}&space;i&space;>&space;0&space;\qquad\qquad\quad\&space;deletion&space;\\&space;&space;&space;s[i,j-1]&space;&plus;&space;\delta(-,w_j)&space;&\text{if&space;}&space;j&space;>&space;0&space;\qquad\qquad\quad\&space;insertion&space;\\&space;&space;s[i-1,j-1]&space;&plus;&space;\delta(v_i,w_j)&space;&\text{if&space;}&space;i&space;>&space;0&space;\text{&space;and&space;}&space;j&space;>&space;0&space;\quad&space;match/mismatch\end{cases}" title="s[i,j] = \max\begin{cases} 0 &\text{if } i=0 \text{ and } j=0 \\ s[i-1,j] + \delta(v_i,-) &\text{if } i > 0 \qquad\qquad\quad\ deletion \\ s[i,j-1] + \delta(-,w_j) &\text{if } j > 0 \qquad\qquad\quad\ insertion \\ s[i-1,j-1] + \delta(v_i,w_j) &\text{if } i > 0 \text{ and } j > 0 \quad match/mismatch\end{cases}" />

Needleman-Wunsch algorithm uses dynamic programming approach converting the alignment problem to a bottom-up filling out dp table problem. 

Time Complexity: <img src="https://latex.codecogs.com/svg.image?O(mn)=O(n^2)" title="O(mn)=O(n^2)" /> 

Space Complexity: <img src="https://latex.codecogs.com/svg.image?O(mn)=O(n^2)" title="O(mn)=O(n^2)" />

<a name="hirschberg"></a>
## Hirschberg
![](Hirschberg_pseudocode.png)

<img src="https://latex.codecogs.com/svg.image?(i,j)" title="(i,j)" /> is the starting point and <img src="https://latex.codecogs.com/svg.image?(i',j')" title="(i',j')" /> is the ending point. 

<img src="https://latex.codecogs.com/svg.image?v" title="v" /> and <img src="https://latex.codecogs.com/svg.image?w" title="w" /> are the input sequences. <img src="https://latex.codecogs.com/svg.image?\delta" title="\delta" /> is the scoring function. 

<img src="https://latex.codecogs.com/svg.image?LAST\_COL" title="LAST\_COL" /> is a modified Needleman-Wunsch algorithm that returns the rightmost column of the dynamic programming table. 

Hirschberg algorithm uses the idea of divide and conquer. For each recursion, we fill out the dynamic programming table <img src="https://latex.codecogs.com/svg.image?(i,j)" title="(i,j)" /> to <img src="https://latex.codecogs.com/svg.image?(i',j')" title="(i',j')" />, find the optimal row in the middle column, divide the problem into two subproblems and conquer them recursively.

Time Complexity: 
<img src="https://latex.codecogs.com/svg.image?O(area&plus;area/2&plus;area/4&plus;\cdots)&space;\le&space;O(2area)&space;=&space;O(2mn)=O(n^2)" title="O(area+area/2+area/4+\cdots) \le O(2area) = O(2mn)=O(n^2)" /> 

Same order but 2 times slower than the Needleman-Wunsch algorithm

Space: <img src="https://latex.codecogs.com/svg.image?O(m)" title="O(m)" /> only store two related columns when calculating the prefix and suffix and only store the prefix column and suffix column when searching for <img src="https://latex.codecogs.com/svg.image?i^*" title="i^*" />.

<a name="analysis"></a>
## Analysis
The codes for analysis and comparison are in *Analysis_Comparison_Hirschberg.ipynb*. To run the analysis, it is better to run the file in Google colab. Otherwise, you might need to install the memory\_profiler library on your own. We use the *profile* function decorator in *memory\_profiler* library, *\%memit* and *\%mprun* in google colab to trace the memory usage and we use the *\%prun* in google colab to trace the running time. The input sequence are generated randomly by *numpy.random.choice* function and the scoring function $\delta$ we use is defined as follows: 

<img src="https://latex.codecogs.com/svg.image?\delta(c,d)&space;=&space;\begin{cases}&space;&space;&space;&space;&space;1&space;&\text{if&space;}&space;c=d\not=-&space;\\&space;&space;&space;&space;&space;-1&space;&\text{if&space;}&space;c\not=d&space;\\&space;&space;&space;&space;&space;\infty&space;&\text{if&space;}&space;c=d=-\end{cases}" title="\delta(c,d) = \begin{cases} 1 &\text{if } c=d\not=- \\ -1 &\text{if } c\not=d \\ \infty &\text{if } c=d=-\end{cases}" />

|                  | Needleman-Wunsch | Hirschberg |
|:----------------:|:----------------:|:----------:|
| Time Complexity  |     <img src="https://latex.codecogs.com/svg.image?O(n^2)" title="O(n^2)" />     |  <img src="https://latex.codecogs.com/svg.image?O(n^2)" title="O(n^2)" />  | 
| Space Complexity |     <img src="https://latex.codecogs.com/svg.image?O(n^2)" title="O(n^2)" />     |   <img src="https://latex.codecogs.com/svg.image?O(n)" title="O(n)" />   |

<a name="run"></a>
## Run
There are two .py files containing the implementation for corresponding algorithms. To run the file, use the following commands in terminal. You can use optional arguments to align the sequences you want. 

    $ python Needleman-Wunsch.py [--help|-h] [--i float] [--d float] [--s float] [--m float] 
                                    [--seqA str] [--seqB str] [--lenA int] [--lenB int]

    $ python Hirschberg.py [--help|-h] [--i float] [--d float] [--s float] [--m float] 
                                    [--seqA str] [--seqB str] [--lenA int] [--lenB int]

    optional arguments:
    -h, --help   show this help message and exit
    --i I        Insertion score
    --d D        Deletion score
    --s S        Substitution score
    --m M        Match score
    --seqA SEQA  The first input sequence
    --seqB SEQB  The second input sequence
    --lenA LENA  The length of the first input sequence
    --lenB LENB  The length of the second input sequence
