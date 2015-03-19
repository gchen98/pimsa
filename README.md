# Introduction #

Add your content here.


# Details #

Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages

# Summary #

PIMSA is an implementation of a reversible jump MCMC algorithm that samples across the space of all possible models in the context of variable selection.  Output from the sampler can be used for inference of a posterior distribution.  Key features of this program are the ability to sample higher-order interactions (the current version is hard-coded at a limit of 5), and incorporate a variety of prior sources of knowledge to inform the variable selection kernel.

# Requirements #

We have tested this software at this moment on a Linux environment only.  Please contact me if you would like help with your environment.

In order for the program to support exploration of a very vast search space, we implemented the program to integrate to a MySQL database.  If you use the yum repository manager this would be a matter of typing (under root):

prompt> yum install mysql-server mysql-devel

Other dependencies are below.  Type at the prompt:

prompt> yum install boost-devel gsl-devel

# Compilation #

Once the above packages are installed, compilation entails entering the directory named src and running make.  You may have to change the first few lines in makefile that declare the variables CFLAGS and LINKFLAGS if your installation of MySQL and GSL are not in standard locations.

# Configuration #

Rather than require all parameters to be specified on a command line, we opted for a XML based configuration file, where concepts can be organized in a more readable fashion via a hierarchy.  The following is the standard template file that is shipped with the program:



&lt;settings&gt;


> 

&lt;database&gt;


> > 

&lt;enable&gt;

true

&lt;/enable&gt;


> > 

&lt;host&gt;

mec.usc.edu

&lt;/host&gt;


> > 

&lt;user&gt;

garyc

&lt;/user&gt;


> > 

&lt;pw&gt;

m3teknik

&lt;/pw&gt;


> > 

&lt;name&gt;

pathway\_mcmc\_folate

&lt;/name&gt;



> 

&lt;/database&gt;


> 

&lt;inputfiles&gt;


> > 

&lt;trait&gt;

phenotypes.txt

&lt;/trait&gt;


> > 

&lt;snplist&gt;

rslist.txt

&lt;/snplist&gt;


> > 

<study\_geno>

genotypes.txt

</study\_geno>


> > 

<study\_env\_var>

envcov.txt

</study\_env\_var>


> > 

&lt;initmodel&gt;

initmodel.txt

&lt;/initmodel&gt;


> > 

<annotation\_files>


> > > 

<z\_matrix>

mygo.z2

</z\_matrix>


> > > 

<a\_matrix>

mygo.a2

</a\_matrix>



> > 

</annotation\_files>


> > 

<endo\_files>


> > > 

<prior\_geno>

genotypes.prior

</prior\_geno>


> > > 

<prior\_env\_var>

envcov.prior

</prior\_env\_var>


> > > 

<prior\_endopheno>

biomarkers.txt

</prior\_endopheno>


> > > 

<selected\_endo>

biomarkers.z

</selected\_endo>



> > 

</endo\_files>



> 

&lt;/inputfiles&gt;


> 

&lt;sampling&gt;


> > 

<use\_endoprior>

true

</use\_endoprior>

 <!--endo,annotation-->
> > 

<use\_a>

true

</use\_a>


> > 

<use\_z>

true

</use\_z>


> > 

&lt;logistic&gt;

true

&lt;/logistic&gt;


> > 

<marginal\_prior>

false

</marginal\_prior>


> > 

<max\_order>

2

</max\_order>


> > 

<main\_effect\_pref>

3

</main\_effect\_pref>


> > 

&lt;iterations&gt;

10

&lt;/iterations&gt;


> > 

<max\_modelsize>

20

</max\_modelsize>



> 

&lt;/sampling&gt;




&lt;/settings&gt;



**Database section**

enable: set to true to use the MySQL instance
host: hostname where the MySQL server resides (highly recommended that it be run locally to minimize latency)
user: the username used to log in
password: the password for the username
name: a database name on the server that the user has privileges to create and drop tables

**Inputfiles**

The formats for these files will be discussed in greater detail in the following section File Format

trait: This file contains either affection status or continuous valued trait data for the study data
snplist: Contains the names of the SNP identifiers in the study data.  The names are not important.  Rather, the length of this file is used to infer how many SNPs to expect in the genotype files
study\_geno: Contains the genotypes across individuals and SNPs in the study.
study\_env\_var: Contains environmental data across individuals for the study.
initmodel: Contains a list of SNPs and/or environmental variables to be used to construct the initial model in the sampler.

Annotation files: This section under the Inputfiles section is relevant when the desired source of priors come from annotation information (e.g. Gene Ontology, SNP mutation class, etc)

z\_matrix: The values in this file are used to directly populate the Z matrix (fixed effects component) in the hierarchical model.
a\_matrix: Correlations between values in this file are used to populate elements in the A matrix (random effects component) of the hierarchical model.

Endophenotype files:  Alternatively, a user may have a second data set of genotypes where an independent set of individuals were typed on the same (or overlapping) markers as in the study.  These individuals also had biomarker measurements made (e.g. gene-expression data).  The files under this section is relevant for such a design.

prior\_geno: This genotype file is analogous to study\_geno but pertains to the individuals in the independent prior dataset.
prior\_env\_var:  Analogous to the study\_env\_var file
prior\_endopheno: Contains a matrix of biomarker values across all individuals and all biomarkers
selected\_endo: Specifies the subset of biomarkers from prior\_endopheno that will be used to generate columns in the Z matrix.

**Sampling section**
The various parameters in this section are used to customize the behavior of the MCMC sampler

use\_endoprior: This true/false flag instructs the program to read in the files appropriate to the type of prior study design desired.
use\_a: if set to false, the random effects component is ignored
use\_z: if set to false, the fixed effects component is ignored
logistic: set this flag to true if the trait file contains binary outcomes (1/0), otherwise set to false (which implies linear regression)
marginal\_prior: If prior odds are to be estimated, one can enable this flag to true, which would then sample from the marginal prior distribution.  This can be used to generate a distribution that is useful in estimating a Bayes Factor from posterior odds and prior odds ratios:
Bayes Factor = [p(H1|data)/p(H0|data)]/[p(H1)/p(H0)]
max\_order: the maximum order of the interactions desired.  For example, 2 denotes only GxG or GxE interactions. 3 would allow 3 way interactions
main\_effect\_pref: This tuning parameter enters the algorithm as the alpha parameter of a Beta (alpha,beta) distribution when determining the probability of adding a main effect to the current model vs an interaction.  Higher numbers give more preference to main effects, where a value of 1 denotes equal probability between choosing a main effect or interaction.
iterations: The total number of MCMC realizations the program should complete before exiting.
max\_modelsize: The kernel that determines whether the next model mutation should be an add or delete considers the current model size and the max model size so that the probability of adding a variable drops sharply as the current model size approaches the max model size value.  Set this value to approximately double of the desired average number of variables in a model.

File formats:

The following briefly describes the file format for each of the files.  Examples are provided in samples/endophenotype and samples/annotation.   Suppose we have N individuals in the study, M markers, E environmental variables, B biomarkers, and P individuals in the prior database.

trait: This file must contain only 1 column and N rows where each value is the affection status or the measured trait of each individual in the study

snplist: This file must contain only 1 column and M rows, where each value in the file can be arbitrary.

study\_geno: This file must contain M columns and N rows.  Note that there are no delimiters in this file.  Each character in the file denotes a genotype for a SNP/individual combination, where the genotypes are coded as 0=missing,1=homozygous wildtype,2=heterozygous,3=homozygous mutant.

study\_env\_var: This optional file contains E columns and N rows, where the columns are delimited by spaces or tabs.  The values are environmental variable data.

initmodel: This file contains 1 column and no greater than M + E rows.  The values are zero-indexed so that zero denotes the first SNP.  To include an environmental covariate, simply add M to the environmental covariate’s index (e.g. the 2nd env. var’s index is (2-1)+M).

z\_matrix: This file must contain M+E rows and any number of tab/space delimited columns.  It is recommended that the number of columns be small for performance reasons.  For example, one column in z\_matrix would generate a Z matrix that includes two columns: the intercept and the values in z\_matrix.

a\_matrix: This file must contain M+E rows and any number of tab/space delimited columns.  The columns are usually coded as 1‘s or 0‘s.  Columns represent concepts.  For example, column 1 might represent pathway 1, column 2 as pathway 2, etc.   A SNP is then connected to another SNP if they both share a 1 in the same column(s).

prior\_geno: This file must contain M columns and P rows.  Note that there are no delimiters in this file.  Each character in the file denotes a genotype for a SNP/individual combination, where the genotypes are coded as 0=missing,1=homozygous wildtype,2=heterozygous,3=homozygous mutant.  IMPORTANT: please note that the order of the columns must match the order specified in study\_geno

prior\_env\_var: This optional file contains E columns and N rows, where the columns are delimited by spaces or tabs.  The values are environmental variable data.  IMPORTANT: please note that the order of the columns must match the order specified in study\_env\_var

prior\_endopheno: This file must contain B tab/space delimited columns and P rows.  The values are endophenotypes measured across all B biomarkers and P individuals.

selected\_endo: This file must contain 1 column and no more than B rows.  It is recommended that the number of columns be small for performance reasons.  For example, one row in selected\_endo would generate a Z matrix that includes two columns: the intercept and the correlations between SNPs (from prior\_geno) and specified biomarkers (from prior\_endopheno).

# Output files #

betas.out: Contains the posterior estimates of beta and its standard errors

models.out: Lists the variables at each sampled model

prior\_betas.out: Contains various bits of information used in calculation of the posterior estimates of beta.

prior\_mean.out: Contains posterior estimates of pi, the fixed effect term in the mixed model defining beta’s distribution
prior\_var.out: Contains the residual variances tau<sup>2 and sigma</sup>2 after accounting for the random effect and fixed effect terms respectively.