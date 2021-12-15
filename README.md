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

<a name="global"></a>
## Global Alignment
Given two strings <img src="https://latex.codecogs.com/svg.image?\boldsymbol{v}&space;\in&space;\Sigma^m" title="\boldsymbol{v} \in \Sigma^m" /> and <img src="https://latex.codecogs.com/svg.image?\boldsymbol{w}&space;\in&space;\Sigma^n" title="\boldsymbol{w} \in \Sigma^n" />, where <img src="https://latex.codecogs.com/svg.image?\Sigma" title="\Sigma" /> is the set of characters, a scoring function <img src="https://latex.codecogs.com/svg.image?\delta:(\Sigma&space;\cup&space;\{-\})\times(\Sigma&space;\cup&space;\{-\})&space;\rightarrow&space;\mathbb{R}" title="\delta:(\Sigma \cup \{-\})\times(\Sigma \cup \{-\}) \rightarrow \mathbb{R}" />, find alignment with maximum score.

An alignment <img src="https://latex.codecogs.com/svg.image?\boldsymbol{A}=[a_{i,j}]" title="\boldsymbol{A}=[a_{i,j}]" /> is a <img src="https://latex.codecogs.com/svg.image?2&space;\times&space;k" title="2 \times k" /> matrix s.t. 
- <img src="https://latex.codecogs.com/svg.image?k&space;=&space;\{max(m,n),\cdots,m&plus;n\}" title="k = \{max(m,n),\cdots,m+n\}" />
- <img src="https://latex.codecogs.com/svg.image?a_{i,j}&space;\in&space;\Sigma&space;\cup&space;\{-\}" title="a_{i,j} \in \Sigma \cup \{-\}" />
- there is no <img src="https://latex.codecogs.com/svg.image?j&space;\in&space;[k]" title="j \in [k]" /> where <img src="https://latex.codecogs.com/svg.image?a_{1,j}=a_{2,j}=-" title="a_{1,j}=a_{2,j}=-" />

<a name="needleman-wunsch"></a>
## Needleman-Wunsch
$$
s[i,j] = max
\begin{cases}
	 0 &\text{if } i=0 \text{ and } j=0 \\
	 s[i-1,j] + \delta(v_i,-) &\text{if } i > 0 &deletion \\ 
	 s[i,j-1] + \delta(-,w_j) &\text{if } j > 0 &insertion \\
	 s[i-1,j-1] + \delta(v_i,w_j) &\text{if } i > 0 \text{ and } j > 0 &match/mismatch
\end{cases}
$$

Needleman-Wunsch algorithm uses dynamic programming approach converting the alignment problem to a bottom-up filling out dp table problem. 

Time Complexity: $O(mn)=O(n^2)$
Space Complexity: $O(mn)=O(n^2)$

<a name="hirschberg"></a>
## Hirschberg
$\textbf{Hirschberg}(i,j,i',j'): \\
    \qquad \text{if } j' - j \gt 1 \\
    \qquad \qquad mid \gets j+\frac{j'-j}{2} \\
    \qquad \qquad prefix \gets \textbf{LAST\_COL}(v[i:i'], w[j:mid], \delta) \\
    \qquad \qquad suffix \gets \textbf{LAST\_COL}(reverse(v[i:i']), reverse(w[mid:j']), \delta) \\
    \qquad \qquad wt = prefix + suffix \\
    \qquad \qquad i^* \gets \underset{i \leq i'' \leq i'}{\arg\max}\ wt[i''] \\
    \qquad \qquad \text{Report} (i^*,mid) \\
    \qquad \qquad \textbf{Hirschberg}(i, j, i^*, mid) \\
    \qquad \qquad \textbf{Hirschberg}(i^*, mid, i', j') \\
$

$(i,j)$ is the starting point and $(i',j')$ is the ending point. 
$v$ and $w$ are the input sequences. $\delta$ is the scoring function.
$LAST\_COL$ is a modified Needleman-Wunsch algorithm that returns the rightmost column of the dynamic programming table. 

Hirschberg algorithm uses the idea of divide and conquer. For each recursion, we fill out the dynamic programming table $(i,j)$ to $(i',j')$, find the optimal row in the middle column, divide the problem into two subproblems and conquer them recursively.

Time Complexity: 
$O(area+area/2+area/4+\cdots) \le O(2area) = O(2mn)=O(n^2)$
Same order but 2 times slower than the Needleman-Wunsch algorithm

Space: $O(m)$ only store two related columns when calculating the prefix and suffix and only store the prefix column and suffix column when searching for $i^*$.

<a name="analysis"></a>
## Analysis
The codes for analysis and comparison are in *Analysis_Comparison_Hirschberg.ipynb*. To run the analysis, it is better to run the file in Google colab. Otherwise, you might need to install the memory\_profiler library on your own. We use the *profile* function decorator in *memory\_profiler* library, *\%memit* and *\%mprun* in google colab to trace the memory usage and we use the *\%prun* in google colab to trace the running time. The input sequence are generated randomly by *numpy.random.choice* function and the scoring function $\delta$ we use is defined as follows: 

$
\delta(c,d) = 
\begin{cases}
     1 &\text{if } c=d\not=- \\
     -1 &\text{if } c\not=d \\
     \infty &\text{if } c=d=-
\end{cases}
$

|                  | Needleman-Wunsch | Hirschberg |
|:----------------:|:----------------:|:----------:|
| Time Complexity  |     $O(n^2)$     |  $O(n^2)$  | 
| Space Complexity |     $O(n^2)$     |   $O(n)$   |

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
