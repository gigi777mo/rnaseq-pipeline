# Design Formulas and Contrasts in DESeq2

## Basic design

```r
~ condition
```

Compares all levels of `condition` against the reference level.

## With batch correction

```r
~ batch + condition
```

This is the most common design when you have known technical batches.

## Interaction (treatment effect depending on genotype)

```r
~ genotype + treatment + genotype:treatment
```

## Setting the reference level

In the sample sheet, make sure the reference group is recognized. In the R script we use:

```r
samples$condition <- relevel(samples$condition, ref = "control")
```

## Multiple contrasts

After running `DESeq(dds)` you can extract different contrasts:

```r
results(dds, contrast = c("condition", "treated", "control"))
results(dds, contrast = c("condition", "treated2", "control"))
```

## Tips

- Always look at the PCA plot colored by both condition **and** batch.
- If batch and condition are completely confounded, you cannot correct for batch.
- For complex designs, consider using `apeglm` or `ashr` for LFC shrinkage.
