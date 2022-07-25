MPN_calc_loglike_for_optim_showrun <- function (params, BioGeoBEARS_run_object, phy, tip_condlikes_of_data_on_each_state, 
          print_optim = TRUE, areas_list = areas_list, states_list = states_list, 
          force_sparse = force_sparse, cluster_already_open = cluster_already_open, 
          return_what = "loglike", calc_ancprobs = FALSE) 
{
  defaults = "\n\tprint_optim=TRUE; areas_list=areas_list; states_list=states_list; force_sparse=force_sparse; cluster_already_open=cluster_already_open; return_what=\"loglike\"; calc_ancprobs=TRUE\n\t"
  if (is.null(states_list) == FALSE) {
    if (is.na(states_list[[1]]) == FALSE) {
      if (states_list[[1]] == "_") {
        states_list[[1]] = NA
      }
    }
  }
  if (is.null(BioGeoBEARS_run_object$printlevel)) {
    BioGeoBEARS_run_object$printlevel = 0
  }
  printlevel = BioGeoBEARS_run_object$printlevel
  traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE
  m = NULL
  jts_matrix = NULL
  BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
  BioGeoBEARS_model_object = params_into_BioGeoBEARS_model_object(BioGeoBEARS_model_object = BioGeoBEARS_model_object, 
                                                                  params = params)
  if (BioGeoBEARS_run_object$rescale_params == TRUE) {
    BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table = BioGeoBEARS_model_object@params_table)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_model_object@params_table
  }
  BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
  d = BioGeoBEARS_model_object@params_table["d", "est"]
  e = BioGeoBEARS_model_object@params_table["e", "est"]
  a = BioGeoBEARS_model_object@params_table["a", "est"]
  b = BioGeoBEARS_model_object@params_table["b", "est"]
  phy$edge.length = phy$edge.length^b
  areas = areas_list
  if ((is.null(BioGeoBEARS_run_object$list_of_distances_mats) == 
       FALSE)) {
    distances_mat = BioGeoBEARS_run_object$list_of_distances_mats[[1]]
  }
  else {
    distances_mat = matrix(1, nrow = length(areas), ncol = length(areas))
  }
  x = BioGeoBEARS_model_object@params_table["x", "est"]
  dispersal_multipliers_matrix = distances_mat^x
  if ((is.null(BioGeoBEARS_run_object$list_of_envdistances_mats) == 
       FALSE)) {
    envdistances_mat = BioGeoBEARS_run_object$list_of_envdistances_mats[[1]]
  }
  else {
    envdistances_mat = matrix(1, nrow = length(areas), ncol = length(areas))
  }
  n = BioGeoBEARS_model_object@params_table["n", "est"]
  dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
    (envdistances_mat^n)
  if ((is.null(BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats) == 
       FALSE)) {
    manual_dispersal_multipliers_matrix = BioGeoBEARS_run_object$list_of_dispersal_multipliers_mats[[1]]
  }
  else {
    manual_dispersal_multipliers_matrix = matrix(1, nrow = length(areas), 
                                                 ncol = length(areas))
  }
  w = BioGeoBEARS_model_object@params_table["w", "est"]
  dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
    (manual_dispersal_multipliers_matrix^w)
  alt_distance_models_abbr = c("WALD", "HNORM")
  tmp_param_names = row.names(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table)
  TF = grepl(pattern = "WALD", x = tmp_param_names)
  if (sum(TF) > 0) {
    tmpname = "WALD_mu"
    TF = grepl(pattern = tmpname, x = tmp_param_names)
    WALD_mu = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, 
                                                                           "est"]
    tmpname = "WALD_lambda"
    TF = grepl(pattern = tmpname, x = tmp_param_names)
    WALD_lambda = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, 
                                                                               "est"]
    tmp_multipliers = statmod::dinvgauss(x = distances_mat, 
                                         mean = WALD_mu, shape = WALD_lambda)
    normalizer = statmod::dinvgauss(x = 0, mean = WALD_mu, 
                                    shape = WALD_lambda)
    tmp_multipliers = tmp_multipliers/normalizer
    dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
      tmp_multipliers
  }
  TF = grepl(pattern = "HNORM", x = tmp_param_names)
  if (sum(TF) > 0) {
    tmpname = "HNORM_theta"
    TF = grepl(pattern = tmpname, x = tmp_param_names)
    HNORM_theta = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[tmpname, 
                                                                               "est"]
    tmp_multipliers = fdrtool::dhalfnorm(x = distances_mat, 
                                         theta = HNORM_theta)
    normalizer = fdrtool::dhalfnorm(x = 0, theta = HNORM_theta)
    tmp_multipliers = tmp_multipliers/normalizer
    dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
      tmp_multipliers
  }
  dmat_times_d = dispersal_multipliers_matrix * matrix(d, 
                                                       nrow = length(areas), ncol = length(areas))
  amat = dispersal_multipliers_matrix * matrix(a, nrow = length(areas), 
                                               ncol = length(areas))
  if ((is.null(BioGeoBEARS_run_object$list_of_area_of_areas) == 
       FALSE)) {
    area_of_areas = BioGeoBEARS_run_object$list_of_area_of_areas[[1]]
  }
  else {
    area_of_areas = rep(1, length(areas))
  }
  u = BioGeoBEARS_model_object@params_table["u", "est"]
  extinction_modifier_list = area_of_areas^(1 * u)
  elist = extinction_modifier_list * rep(e, length(areas))
  if (is.null(BioGeoBEARS_run_object$custom_Qmat_fn_text) == 
      TRUE) {
    if (traitTF == FALSE) {
      Qmat = rcpp_states_list_to_DEmat(areas_list = areas_list, 
                                       states_list = states_list, dmat = dmat_times_d, 
                                       elist = elist, amat = amat, include_null_range = BioGeoBEARS_run_object$include_null_range, 
                                       normalize_TF = TRUE, makeCOO_TF = force_sparse)
    }
    if (traitTF == TRUE) {
      numstates_geogtrait = ncol(tip_condlikes_of_data_on_each_state)
      res = modify_Qmat_with_trait(Qmat = NULL, BioGeoBEARS_run_object, 
                                   numstates_geogtrait = numstates_geogtrait, areas_list = areas_list, 
                                   states_list = states_list, dispersal_multipliers_matrix = dispersal_multipliers_matrix, 
                                   elist = elist, force_sparse = force_sparse)
      Qmat = res$Qmat
      m = res$m
      if (is.null(BioGeoBEARS_run_object$jts_txt_matrix) == 
          FALSE) {
        jts_txt_matrix = BioGeoBEARS_run_object$jts_txt_matrix
        jts_matrix = matrix(data = 0, nrow = nrow(jts_txt_matrix), 
                            ncol = ncol(jts_txt_matrix))
        TF_matrix = matrix(data = TRUE, nrow = nrow(jts_txt_matrix), 
                           ncol = ncol(jts_txt_matrix))
        diag(TF_matrix) = FALSE
        jts_txt_params = c(jts_txt_matrix[TF_matrix])
        jts_txt_params
        for (jts_i in 1:nrow(jts_txt_matrix)) {
          diag_val = 1
          for (jts_j in 1:ncol(jts_txt_matrix)) {
            if (jts_i == jts_j) {
              (next)()
            }
            jts_txt = jts_txt_matrix[jts_i, jts_j]
            newval = as.numeric(BioGeoBEARS_model_object@params_table[jts_txt, 
                                                                      "est"])
            jts_matrix[jts_i, jts_j] = newval
            diag_val = 1 - newval
          }
          jts_matrix[jts_i, jts_i] = diag_val
        }
      }
    }
  }
  else {
    cat("\n\nNOTE: BioGeoBEARS is using a custom Qmat-generating function.\n\n")
    eval(parse(text = BioGeoBEARS_run_object$custom_Qmat_fn_text))
  }
  spPmat_inputs = get_spPmat_inputs_from_BGB(BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                             states_list = states_list, dispersal_multipliers_matrix = dispersal_multipliers_matrix)
  if (is.null(BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state) == 
      TRUE) {
    if (BioGeoBEARS_run_object$use_detection_model == TRUE) {
      numareas = length(areas)
      detects_df = BioGeoBEARS_run_object$detects_df
      controls_df = BioGeoBEARS_run_object$controls_df
      mean_frequency = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mf", 
                                                                                    "init"]
      dp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["dp", 
                                                                        "init"]
      fdp = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["fdp", 
                                                                         "init"]
      tip_condlikes_of_data_on_each_state = tiplikes_wDetectionModel(states_list_0based_index = states_list, 
                                                                     phy = phy, numareas = numareas, detects_df = detects_df, 
                                                                     controls_df = controls_df, mean_frequency = mean_frequency, 
                                                                     dp = dp, fdp = fdp, null_range_gets_0_like = TRUE, 
                                                                     return_LnLs = TRUE, relative_LnLs = TRUE, exp_LnLs = TRUE, 
                                                                     error_check = TRUE)
      if (is.null(BioGeoBEARS_run_object$prior_by_range_size) == 
          FALSE) {
        for (iii in 1:nrow(tip_condlikes_of_data_on_each_state)) {
          tip_condlikes_of_data_on_each_state[iii, ] = tip_condlikes_of_data_on_each_state[iii, 
          ] * BioGeoBEARS_run_object$prior_by_range_size
        }
      }
    }
  }
  else {
    tip_condlikes_of_data_on_each_state = BioGeoBEARS_run_object$tip_condlikes_of_data_on_each_state
  }
  if (print_optim == TRUE) {
  }
  if (calc_ancprobs == FALSE) {
    fixnode = BioGeoBEARS_run_object$fixnode
    fixlikes = BioGeoBEARS_run_object$fixlikes
    ttl_loglike = calc_loglike_sp(tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                  phy = phy, Qmat = Qmat, spPmat = NULL, return_what = "loglike", 
                                  probs_of_states_at_root = NULL, sparse = force_sparse, 
                                  use_cpp = TRUE, input_is_COO = force_sparse, spPmat_inputs = spPmat_inputs, 
                                  cppSpMethod = 3, printlevel = BioGeoBEARS_run_object$printlevel, 
                                  cluster_already_open = cluster_already_open, calc_ancprobs = FALSE, 
                                  fixnode = fixnode, fixlikes = fixlikes, include_null_range = BioGeoBEARS_run_object$include_null_range, 
                                  states_allowed_TF = NULL, m = m, jts_matrix = jts_matrix, 
                                  BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                  on_NaN_error = BioGeoBEARS_run_object$on_NaN_error)
    ttl_loglike
    if (print_optim == TRUE) {
      LnL = ttl_loglike
      outvars = adf(t(c(BioGeoBEARS_model_object@params_table$est, 
                        LnL)))
      names(outvars) = c(rownames(BioGeoBEARS_model_object@params_table), 
                         "LnL")
      print(round(outvars, 3))
    }
    return(ttl_loglike)
  }
  else {
    fixnode = BioGeoBEARS_run_object$fixnode
    fixlikes = BioGeoBEARS_run_object$fixlikes
    model_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state = tip_condlikes_of_data_on_each_state, 
                                    phy = phy, Qmat = Qmat, spPmat = NULL, return_what = "all", 
                                    probs_of_states_at_root = NULL, sparse = force_sparse, 
                                    use_cpp = TRUE, input_is_COO = force_sparse, spPmat_inputs = spPmat_inputs, 
                                    printlevel = BioGeoBEARS_run_object$printlevel, 
                                    cluster_already_open = cluster_already_open, calc_ancprobs = TRUE, 
                                    fixnode = fixnode, fixlikes = fixlikes, include_null_range = BioGeoBEARS_run_object$include_null_range, 
                                    states_allowed_TF = NULL, m = m, jts_matrix = jts_matrix, 
                                    BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object, 
                                    on_NaN_error = BioGeoBEARS_run_object$on_NaN_error)
    return(model_results)
  }
  print(paste0("Currently running model ", resfnPLOT, " in ", trfnPH[[1]], " phylogeny.", sep = ""))
}
