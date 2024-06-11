## ----model.matrix---------------------------------------------------------------------
## ?model.matrix
mat <- with(trees, model.matrix(log(Volume) ~ log(Height) + log(Girth)))
mat
colnames(mat)


## ----lm_example-----------------------------------------------------------------------
summary(lm(log(Volume) ~ log(Height) + log(Girth), data = trees))


## ----EMM_example1---------------------------------------------------------------------
## Load ExploreModelMatrix
library("ExploreModelMatrix")

## Example data
(sampleData <- data.frame(
    genotype = rep(c("A", "B"), each = 4),
    treatment = rep(c("ctrl", "trt"), 4)
))

## Let's make the visual aids provided by ExploreModelMatrix
vd <- ExploreModelMatrix::VisualizeDesign(
    sampleData = sampleData,
    designFormula = ~ genotype + treatment,
    textSizeFitted = 4
)

## Now lets plot these images
cowplot::plot_grid(plotlist = vd$plotlist)


## ----EMM_example1_interactive, eval = FALSE-------------------------------------------
## ## We are using shiny again here
app <- ExploreModelMatrix(
    sampleData = sampleData,
    designFormula = ~ genotype + treatment
 )
if (interactive()) shiny::runApp(app)




########### additional codes in the lecture from Leo ##########
devel > df <- data.frame(genomic_background = rep(c("WT", "mutA"), 2))
devel > df
genomic_background
1                 WT
2               mutA
3                 WT
4               mutA
devel > model.matrix(~ genomic_background, data = df)
(Intercept) genomic_backgroundWT
1           1                    1
2           1                    0
3           1                    1
4           1                    0
attr(,"assign")
[1] 0 1
attr(,"contrasts")
attr(,"contrasts")$genomic_background
[1] "contr.treatment"

devel > class(df$genomic_background)
[1] "character"
devel >
    devel >
    devel > factor(df$genomic_background)
[1] WT   mutA WT   mutA
Levels: mutA WT
devel > unique(df$genomic_background)
[1] "WT"   "mutA"
devel > sort(unique(df$genomic_background))
[1] "mutA" "WT"
devel > sort(unique(df$genomic_background))[1]
[1] "mutA"
devel > factor(df$genomic_background, levels = c("WT", "mutA"))
[1] WT   mutA WT   mutA
Levels: WT mutA
devel > df$genomic_factor <- factor(df$genomic_background, levels = c("WT", "mutA"))
devel > model.matrix(~ genomic_factor, data = df)
(Intercept) genomic_factormutA
1           1                  0
2           1                  1
3           1                  0
4           1                  1
attr(,"assign")
[1] 0 1
attr(,"contrasts")
attr(,"contrasts")$genomic_factor
[1] "contr.treatment"

devel >
    devel > df$treatment <- rep(c("treatment", "notreatment"), 2)
devel > df
genomic_background genomic_factor   treatment
1                 WT             WT   treatment
2               mutA           mutA notreatment
3                 WT             WT   treatment
4               mutA           mutA notreatment
devel > model.matrix(~ treatment, data = df)
(Intercept) treatmenttreatment
1           1                  1
2           1                  0
3           1                  1
4           1                  0
attr(,"assign")
[1] 0 1
attr(,"contrasts")
attr(,"contrasts")$treatment
[1] "contr.treatment"


####################################


























