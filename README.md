# LMSstat
Automation of statistical test with an identical data input aiming to reduce arduous work searching for packages and changing data input.

The package includes

* Simple Statistics :u-test, t-test, post-hocs of Anova and Kruskal-Wallis with FDR adjusted values

* Bar, Box, Dot, Violin plots with significance (u-test, t-test, post-hocs of Anova and Kruskal-Wallis)

* Scaling & Transformation

* Normality check (Shapiro Wilk test)

* Scheirer–Ray–Hare Test

* Volcano plot

* Heatmap

* PERMANOVA

* NMDS

* PCA

* PCoA

## Instructions

### Installation

#### Download R
https://cran.r-project.org/bin/windows/base/
#### Download R Studio
https://www.rstudio.com/products/rstudio/download/
#### Download Rtools
https://cran.r-project.org/bin/windows/Rtools/

#### Download package in R

```
install.packages("devtools")

devtools::install_github("SNUFML/LMSstat")

# devtools::install_github("SNUFML/LMSstat", ref = "develop")

library(LMSstat)
```


### Basic structure of the Data
#### Used in

* Simple statistics
* Barplot, Boxplot, Dotplot
* Volcano plot
* Scheirer–Ray–Hare Test
* PERMANOVA
* NMDS
* PCA
* Scaling & Transformation
* Normality check (Shapiro Wilk test)
* Heatmap


```
#Sample Data provided within the package

data("Data")

# Uploading your own Data

setwd("~")

Data <- read.csv("statT.csv", check.names = FALSE)

# When uploading an Excel file
Data <- as.data.frame(readxl::read_excel("statT.xlsx", sheet = 2))

# When your data is large
Data <- as.data.frame(data.table::fread("statT.csv", check.names = FALSE, header = TRUE, integer64 = "double"))
```

## Datafile needs to follow the following format

<p align="center">  
<img src="https://user-images.githubusercontent.com/87584921/183798448-03f22df1-aa29-4a38-89b4-014be7a9078f.png" width="750" height="400" alt="statT.csv">
</p>
<p align="center">statT.csv</p>


#### Used in {id="used-in_1"}

* PERMANOVA

```
# Sample Data provided within the package
data("Classification")

# Uploading your own Data
Classification <- read.csv("statT_G.csv", header = F)
```

<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125877154-d01ad8b0-25cd-448b-905d-749a6cc93552.PNG" width="500" height="400" alt="statT_G.csv">
</p>
<p align="center">statT_G.csv</p>


### Univariate statistics

```
# Use `parallel = TRUE` for datasets with more than 1000 metabolites.
Statfile <- All_stats(Data,Adjust_p_value = TRUE, Adjust_method = "BH", parallel = FALSE) 

Statfile <- Allstats(Data,Adjust_p_value = TRUE, Adjust_method = "BH") # Previous version using for-loop
```
##### Adjustable parameters {id="adjustable-parameters_1"}

* Adjust_p_value = T # Set True if adjustment is needed

* Adjust_method = F # Adjustment methods frequently used. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

```
head(Statfile[["Result"]]) # includes all statistical results

write.csv(Statfile[["Result"]], "p_value_result.csv")  # Write csv with all the p-value included

install.packages("writexl")

writexl::write_xlsx(list(mysheet = Data), "p_value_result.xlsx")  # Write xlsx with all the p-value included

data.table::fwrite(as.data.frame(Statfile[["Result"]]), "p_value_result.csv", row.names = TRUE)
```
### Plots - DO NOT USE

```
# Makes a subdirectory and saves box plots for all the variables
# Boxplot(Statfile,asterisk = "u_test") 

# Makes a subdirectory and saves dot plots for all the variables
# Dotplot(Statfile,asterisk = "t_test") 

# Makes a subdirectory and saves bar plots for all the variables
# Barplot(Statfile,asterisk = "Scheffe")

# Makes a subdirectory and saves violin plots for all the variables
# Violinplot(Statfile,asterisk = "Scheffe")
```


&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**Boxplot(Statfile)**&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**Dotplot(Statfile)**
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125714687-1908d6eb-b2bd-4e25-8ef0-62c24466c32a.png" width="350" height="250" alt="Boxplot">
<img src="https://user-images.githubusercontent.com/77651662/125714704-d7dab67e-03c0-4e35-b86a-36723f7c63de.png" width="350" height="250" alt="Dotplot">
</p>

&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**Barplot(Statfile)**&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**Violinplot(Statfile)**
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125715925-0878ec77-30bf-4859-8e56-316d98b6d520.jpg" width="350" height="250" alt="Barplot">
<img src="https://user-images.githubusercontent.com/77651662/141611647-1bc0f049-fcd1-4cb4-a713-5aa30f24dcb1.png" width="350" height="250" alt="Violinplot">
</p>

##### Adjustable parameters {id="adjustable-parameters_2"}

* asterisk = "t_test"   #c("Dunn","Scheffe","u_test","t_test")
* significant_variable_only = F  # If set to TRUE, insignificant results will not be plotted
* color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099") # Colors for the plots
* legend_position = "none" #  "none","left","right","bottom","top"
* order = NULL # Order of the groups c("LAC","LUE","WEI","SDF","HGH","ASH")
* tip_length = 0.01 # significance tip length
* label_size = 2.88 # significance label size
* step_increase = 0.05 #significance step increase
* width = 0.3 # box width ; size = 3 # dot size
* fig_width = NA #figure size
* fig_height = NA #figure size
* Y_text = 12 # Y title size
* X_text = 10 # X text size
* Y_lab = 10 #y axis text size
* T_size = 15 # Title size
* sig_int = c(0.1,0.05) # significance interval

### Scaling & Transformation

```
scaled_data <- D_tran(Data, param = "Auto")
```
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**Raw_Data**&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**Scaled_Data**
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/126724661-1e2ee121-ecfe-41eb-a52a-e3dca128c12c.PNG" width="400" height="200" alt="scaled_original">
<img src="https://user-images.githubusercontent.com/77651662/126724676-44c4eaac-c007-4eaf-9962-0b739291adb6.PNG" width="400" height="200" alt="scaled">
</p>

##### Adjustable parameters {id="adjustable-parameters_3"}

* param = "None" # "None","Auto","log10","Pareto"

* save = F  #Set true if datafile is to be saved

### Normality check

```
#Shapiro Wilk test

Result <- Norm_test(Data)

write.csv(Result, "Normality_test_Result.csv")
```
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/126578553-f2d583db-d649-4db7-92bb-8cfa9433904a.PNG" width="400" height="350" alt="normality test">
</p>

### Scheirer–Ray–Hare Test - do not use

```
# csv files including significant variables (Multilevel, Group, interaction) and a Venn diagram are downloaded
SRH(Data)
```
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/130378279-fb49da48-dc7b-418e-88e7-db046364a4cc.PNG" width="400" height="350" alt="SRH">
</p>

##### Adjustable parameters {id="adjustable-parameters_4"}

* Adjust_p_value = T # Set True if adjustment is needed
* Adjust_method = "BH" # Adjustment methods frequently used. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")


### Volcano plot

```
# Makes a subdirectory and saves Volcano plots for different combination of groups
Test <- All_stats(Data)
Volcano(Test, asterisk = "t-test")
```
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/127811109-0cb273ba-3c01-4679-9cb3-c510f3176592.png" width="400" height="400" alt="volcano">
<img src="https://user-images.githubusercontent.com/77651662/127816379-409d8630-eaf2-45aa-8d77-7540e5542ebd.png" width="400" height="400" alt="volcano_named">
</p>


##### Adjustable parameters {id="adjustable-parameters_5"}

* asterisk = "t-test" #statistics inheriting from Allstats "Scheffe", "t-test", "u-test", "Dunn"
* reverse = T # T, F reverse the direction of fold change
* fig_width = NA #figure size
* fig_height = NA #figure size
* FC_log = 2 # Fold change log transformation value
* pval_log = 10 #p_value log transformation value
* dotsize = 3 #dotsize
* x_limit = c(-2,2) #x axis limt
* y_limit =c(0,6) #y axis limit
* pval_intercept = 0.05 # intercept for identification
* sig_label = T # T,F label significant variables
* color=c("#FF3300","#FF6600","#FFCC00") #colors used for ggplots.
* fixed_limit = F #whether the limit should be fixed or not T, F
* max_overlap = 20 #maximum overlap for labels
* FC_range = c(-1.5,1.5) #significant fold change range

### Heatmap

```
# Makes a subdirectory and saves Heatmap

scaled_data <- D_tran(Data, param = "Auto")

Heatmap(scaled_data) #data inheriting from D_tran
```
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/126421942-247031e1-9f90-452f-b5bd-fdf5bdf5c058.PNG" width="750" height="400" alt="Heatmap">
</p>


##### Adjustable parameters {id="adjustable-parameters_6"}
* col =c("green", "white", "red") # colors for heatmap
* col_lim = c(-3, 0, 3) # color boundaries
* reverse = T # T,F Reverse column and rows
* distance = "pearson" # Distance matrix for HCA "pearson", "manhattan","euclidean","spearman","kendall" ,
* rownames = T # T,F
* colnames = T # T,F
* Hsize = (3,6) # Width & Height c(a,b)
* g_legend = "Group" # Annotation legend title
* h_legend = "Color Key" # Heatmap legend title
* Title ="Title" # Title
* T_size = 10 # Title text size
* R_size = 3 # row text size
* C_size = 3 # column text size
* Gcol =c("ASD" = "black","HGH"="red","LAC"="blue","LUE" ="grey","SDF" = "yellow","WEI"="green") # Color for top_annotation bar
* dend_h = 0.5 #dendrite height
* a_h = 0.2 # top annotation hegiht

### Multivariate statistics
#### PERMANOVA
```
data("Data")

data("Classification") 
```

#### Single factor
PERMANOVA done with the Group column
```
Indiv_Perm(Data) # The group information is treated as a factor
```

<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125876315-c51448c9-eef1-4299-b85e-34e62a99bea8.PNG" width="750" height="225" alt="PERMANOVA-single">
</p>

#### Multiple Factors

Loops PERMANOVA over different classes provided by Classification

```
Result <- Multi_Perm(Data, Classification) # The group information is treated as factors
```
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125875797-89b696dd-f2a7-4ff3-aec9-74ac4b9075f9.PNG" width="750" height="400" alt="PERMANOVA-multi">
</p>



##### Adjustable parameters {id="adjustable-parameters_7"}
* method = Dissimilarity index c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord")

#### NMDS
```
# Makes a subdirectory and saves NMDS plots for all of the distance metrics
NMDS(Data, methods = c("manhattan", "bray", "euclidean"))
```
<p align="center">
<img src=https://user-images.githubusercontent.com/77651662/125900616-c0c6728d-0b3a-445b-bf41-e32be766924f.png width="600" height="500" alt="NMDS">
</p>
<p align="center">NMDS plot with bray distance and p-value from PERMANOVA</p>

##### Adjustable parameters {id="adjustable-parameters_8"}
* methods = Dissimilarity index c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord")

* color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099") # Colors for the plots
* legend_position = "none" #  "none","left","right","bottom","top"
* fig_width = NA #figure size
* fig_height = NA #figure size
* names = F # used to indicate sample names
* dotsize = 3 # dotsize
* labsize = 3 # label size
#### PCA
```
# Makes a subdirectory and saves PCA plot
PCA(Data, components = c(1, 2), legend_position = "none")
```
<p align="center">
<img src=https://user-images.githubusercontent.com/77651662/126108588-cfe688c6-2c90-485d-a197-ef1eb7a82cb5.png width="600" height="500" alt="PCA">
</p>
<p align="center">PCA plot with selected components</p>

##### Adjustable parameters {id="adjustable-parameters_9"}
* color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099") # Colors for the plots
* legend_position = "none" #  "none","left","right","bottom","top"
* fig_width = NA #figure size
* fig_height = NA #figure size
* components = c(1,2) # selected components
* names = F # used to indicate sample names
* dotsize = 3 # dotsize
* labsize = 3 # label size
* ellipse = T # T  or F to show ellipse
#### PCoA
```
# Makes a subdirectory and saves PCoA plot
PCoA(Data, components = c(1, 2), methods = c("bray", "manhattan"))
```
<p align="center">
<img src=https://user-images.githubusercontent.com/77651662/127433788-7aa75a05-3559-4bd1-9504-c1234c5905d4.png width="600" height="500" alt="PCoA">
</p>
<p align="center">PCoA plot with selected components</p>

##### Adjustable parameters {id="adjustable-parameters_10"}
* color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099") # Colors for the plots
* legend_position = "none" #  "none","left","right","bottom","top"
* fig_width = NA #figure size
* fig_height = NA #figure size
* components = c(1,2) # selected components
* names = F # used to indicate sample names
* dotsize = 3 # dotsize
* labsize = 3 # label size
* ellipse = T # T  or F to show ellipse
* methods = Dissimilarity index c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq",chord")
