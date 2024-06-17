# agriutilities (development version)

  * `extract_rcov()` does not require any other parameter but the model.
  * `ic_reml_asr()` implementation of the information criteria proposed by Verbyla (2019).
  * `ic_reml_spt()` implementation of the information criteria proposed by Verbyla (2019) for SpATS models.
  * Add `h_cullis_spt()` to calculate the generalized heritability proposed by Cullis (2006).
  * `remove_outliers` FALSE by default in `single_trial_analysis()`.
  * Minor changes in `check_design_met()`.
  * Fix absolute value of the mean in CV for `single_trial_analysis()`.

# agriutilities 1.2.0

  * Add `extract_rcov()` to extract variance-covariance matrices.
  * Add `method` to be passed to cor function in `gg_cor()`.
  * Change style for naming the parameters in `gg_cor()`.
  * Minor changes in the writing style.
  * Add `reorder` a logical value to reorder by a Hierarchical Clustering in `covcor_heat()`.

# agriutilities 1.1.0

## Minor changes/fixes

  * Add `type = "spatial"` to `plot.smaAgri()` method.
  * Add restriction when identifying row-column designs.
  * Add minimum and maximum value when returning summary in `check_design_met()`.
  * Add S3 method to plot an object of class metAgri `plot.metAgri()`.
  * Add S3 method to plot an object of class smaAgri `plot.smaAgri()`.
  * Fix argument not working in `covcor_heat()`.

# agriutilities 1.0.0

