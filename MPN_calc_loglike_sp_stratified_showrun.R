MPN_calc_loglike_sp_stratified_showrun <- function (tip_condlikes_of_data_on_each_state, phy, Qmat = NULL, 
          spPmat = NULL, min_branchlength = 1e-06, return_what = "loglike", 
          probs_of_states_at_root = NULL, rootedge = TRUE, sparse = FALSE, 
          printlevel = 0, use_cpp = TRUE, input_is_COO = FALSE, spPmat_inputs = NULL, 
          cppSpMethod = 3, cluster_already_open = NULL, calc_ancprobs = FALSE, 
          include_null_range = TRUE, fixnode = NULL, fixlikes = NULL, 
          inputs = inputs, allareas = allareas, all_states_list = all_states_list, 
          return_condlikes_table = FALSE, calc_TTL_loglike_from_condlikes_table = TRUE) 
{
  defaults = "\n\tQmat=NULL; spPmat=NULL; min_branchlength=0.000001; return_what=\"loglike\"; probs_of_states_at_root=NULL; rootedge=FALSE; sparse=FALSE; printlevel=1; use_cpp=TRUE; input_is_COO=FALSE; spPmat_inputs=NULL; cppSpMethod=3; cluster_already_open=NULL; calc_ancprobs=FALSE; include_null_range=TRUE; fixnode=NULL; fixlikes=NULL; inputs=inputs; allareas=allareas; all_states_list=all_states_list; return_condlikes_table=FALSE; calc_TTL_loglike_from_condlikes_table=TRUE\n\t"
  defaults = "\n\tmaxareas = 4\n\tinclude_null_range = TRUE\n\tphy = read.tree(inputs$trfn)\n\ttipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(inputs$geogfn))\n\ttip_condlikes_of_data_on_each_state = tipranges_to_tip_condlikes_of_data_on_each_state(tipranges, phy, maxareas=maxareas, include_null_range=include_null_range)\n\t\n\tallareas = getareas_from_tipranges_object(tipranges)\n\tall_states_list = rcpp_areas_list_to_states_list(areas=allareas, include_null_range=TRUE, maxareas=maxareas)\n\t\n\ttmpres = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=0.000001, return_what=\"all\", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=FALSE, include_null_range=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=FALSE)\n\ttmpres\n\t\n\tmin_branchlength=0.000001\n\tinclude_null_range=TRUE\n\tprintlevel=0\n\tcppSpMethod=3\n\treturn_condlikes_table=TRUE\n\tcalc_TTL_loglike_from_condlikes_table=TRUE\n\tcalc_ancprobs=TRUE\n"
  defaults = "\n\ttmpres = calc_loglike_sp_stratified(tip_condlikes_of_data_on_each_state, phy, Qmat=NULL, spPmat=NULL, min_branchlength=0.000001, return_what=\"all\", probs_of_states_at_root=NULL, rootedge=TRUE, sparse=FALSE, printlevel=0, use_cpp=TRUE, input_is_COO=FALSE, spPmat_inputs=NULL, cppSpMethod=3, cluster_already_open=NULL, calc_ancprobs=TRUE, include_null_range=TRUE, fixnode=NULL, fixlikes=NULL, inputs=inputs, allareas=allareas, all_states_list=all_states_list, return_condlikes_table=TRUE, calc_TTL_loglike_from_condlikes_table=TRUE)\n"
  BioGeoBEARS_run_object = inputs
  if (is.null(inputs$printlevel)) {
    inputs$printlevel = 0
  }
  printlevel = inputs$printlevel
  traitTF = is.null(BioGeoBEARS_run_object$trait) == FALSE
  if (traitTF == TRUE) {
    trait_Pmat_txt = BioGeoBEARS_run_object$trait_Pmat_txt
    num_trait_states = ncol(trait_Pmat_txt)
  }
  m = NULL
  jts_matrix = NULL
  BioGeoBEARS_model_object = BioGeoBEARS_run_object$BioGeoBEARS_model_object
  if (BioGeoBEARS_run_object$rescale_params == TRUE) {
    BioGeoBEARS_model_object@params_table = unscale_BGB_params(scaled_params_table = BioGeoBEARS_model_object@params_table)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = BioGeoBEARS_model_object@params_table
  }
  BioGeoBEARS_model_object = calc_linked_params_BioGeoBEARS_model_object(BioGeoBEARS_model_object)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object = BioGeoBEARS_model_object
  inputs$BioGeoBEARS_model_object = BioGeoBEARS_model_object
  if (!is.null(fixnode)) {
    if ((is.null(dim(fixlikes)) == TRUE) && (length(fixnode) == 
                                             1)) {
      pass_fixlikes = TRUE
    }
    else {
      if ((dim(fixlikes)[1]) == length(fixnode)) {
        pass_fixlikes = TRUE
        if ((all(c(order(fixnode) == 1:length(fixnode)))) == 
            TRUE) {
          pass_fixlikes = TRUE
        }
        else {
          pass_fixlikes = FALSE
          error_msg = "ERROR in calc_loglike_sp_stratified(): \n             Multiple nodes in 'fixnode' MUST be sorted in increasing order.\n"
          cat(error_msg)
          stop(error_msg)
        }
      }
      else {
        pass_fixlikes = FALSE
        error_msg = "ERROR in calc_loglike_sp_stratified(): Either:\n             (1) fixnode must be a single node number, and fixlikes must be a vector, or\n             (2) fixlikes like must be a matrix with the # of rows equal to length(fixnode).\n"
        cat(error_msg)
        stop(error_msg)
      }
    }
  }
  if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == 
                                           TRUE)) {
    names_in_inputs = names(inputs)
    if (("master_table" %in% names_in_inputs) == TRUE) {
      condlikes_table = matrix(data = 0, nrow = nrow(inputs$master_table), 
                               ncol = ncol(tip_condlikes_of_data_on_each_state))
      tmprownums = nrow(tip_condlikes_of_data_on_each_state)
      condlikes_table[1:tmprownums, ] = tip_condlikes_of_data_on_each_state
      if (calc_ancprobs == TRUE) {
        if (traitTF == FALSE) {
          relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = matrix(data = 0, 
                                                                                           nrow = nrow(inputs$master_table), ncol = length(all_states_list))
        }
        else {
          relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = matrix(data = 0, 
                                                                                           nrow = nrow(inputs$master_table), ncol = length(all_states_list) * 
                                                                                             num_trait_states)
        }
      }
    }
    else {
      cat("\n\nWARNING: in 'calc_loglike_sp_stratified()', you set 'return_condlikes_table=TRUE'\n\n\t\t\tand/or calc_TTL_loglike_from_condlikes_table=TRUE, but this requires that\n\n\t\t\t'inputs$master_table' be available from the 'section_the_tree()' function. Try\n\t\t\t\nrunning 'inputs=section_the_tree(inputs, make_master_table=TRUE).\n", 
          sep = "")
      cat("\nAs a result, we are setting return_condlikes_table=FALSE\n\n", 
          sep = "")
      return_condlikes_table = FALSE
    }
  }
  if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 
      1) {
    num_timeperiods = 1
  }
  else {
    timeperiods = inputs$timeperiods
    num_timeperiods = length(timeperiods)
  }
  allareas = allareas
  allareas_list = seq(0, length(allareas) - 1, 1)
  areas = allareas_list
  all_states_list = all_states_list
  BioGeoBEARS_model_object = inputs$BioGeoBEARS_model_object
  force_sparse = sparse
  tip_relative_probs_of_each_state = tip_condlikes_of_data_on_each_state/rowSums(tip_condlikes_of_data_on_each_state)
  tip_relative_probs_of_each_state
  current_tip_relative_probs_of_each_state = tip_relative_probs_of_each_state
  current_condlikes_row = nrow(current_tip_relative_probs_of_each_state)
  numnodes = phy$Nnode + length(phy$tip.label)
  all_relative_probs_of_each_state = matrix(0, ncol = ncol(tip_condlikes_of_data_on_each_state), 
                                            nrow = (numnodes * length(timeperiods)))
  all_condlikes_of_each_state = matrix(0, ncol = ncol(tip_condlikes_of_data_on_each_state), 
                                       nrow = (numnodes * length(timeperiods)))
  all_relative_probs_of_each_state[1:current_condlikes_row, 
  ] = current_tip_relative_probs_of_each_state
  all_condlikes_of_each_state[1:current_condlikes_row, ] = current_tip_relative_probs_of_each_state
  previous_timepoint = 0
  original_phy = phy
  phy_as_it_is_chopped_down = original_phy
  for (i in 1:num_timeperiods) {
    d = BioGeoBEARS_model_object@params_table["d", "est"]
    e = BioGeoBEARS_model_object@params_table["e", "est"]
    a = BioGeoBEARS_model_object@params_table["a", "est"]
    user_specified_constraints_on_states_list_TF = FALSE
    states_allowed_TF1 = rep(TRUE, length(all_states_list))
    states_allowed_TF2 = rep(TRUE, length(all_states_list))
    states_allowed_TF3 = rep(TRUE, length(all_states_list))
    if ((is.null(inputs$list_of_areas_allowed_mats) == FALSE)) {
      user_specified_constraints_on_states_list_TF = TRUE
    }
    if ((is.null(inputs$list_of_areas_adjacency_mats) == 
         FALSE)) {
      user_specified_constraints_on_states_list_TF = TRUE
    }
    if ((is.null(inputs$lists_of_states_lists_0based) == 
         FALSE)) {
      user_specified_constraints_on_states_list_TF = TRUE
    }
    if (user_specified_constraints_on_states_list_TF == 
        TRUE) {
      if (is.null(inputs$lists_of_states_lists_0based) == 
          TRUE) {
        errortxt = paste0("STOP ERROR in calc_loglike_sp_stratified(): User has specified areas_allowed or area_adjacency constraints, but 'lists_of_states_lists_0based' has not been added to the BioGeoBEARS_run_object.")
        cat("\n\n")
        cat(errortxt)
        cat("\n\n")
        stop(errortxt)
      }
      if ((is.null(inputs$list_of_areas_allowed_mats) == 
           FALSE)) {
        areas_allowed_mat = inputs$list_of_areas_allowed_mats[[i]]
        states_allowed_TF1 = sapply(X = all_states_list, 
                                    FUN = check_if_state_is_allowed, areas_allowed_mat)
        if (include_null_range == TRUE) {
          states_allowed_TF1[1] = TRUE
        }
      }
      if ((is.null(inputs$list_of_areas_adjacency_mats) == 
           FALSE)) {
        areas_adjacency_mat = inputs$list_of_areas_adjacency_mats[[i]]
        states_allowed_TF2 = sapply(X = all_states_list, 
                                    FUN = check_if_state_is_allowed_by_adjacency, 
                                    areas_adjacency_mat)
        if (include_null_range == TRUE) {
          states_allowed_TF2[1] = TRUE
        }
      }
      if ((is.null(inputs$lists_of_states_lists_0based) == 
           FALSE)) {
        states_allowed_TF3 = all_states_list %in% inputs$lists_of_states_lists_0based[[i]]
        if (include_null_range == TRUE) {
          states_allowed_TF3[1] = TRUE
        }
      }
      states_allowed_TF = ((states_allowed_TF1 + states_allowed_TF2 + 
                              states_allowed_TF3) == 3)
      inputs$lists_of_states_lists_0based[[i]] = all_states_list[states_allowed_TF]
    }
    else {
      pass = 1
      states_allowed_TF = rep(TRUE, length(all_states_list))
    }
    states_to_use_TF = rep(TRUE, length(all_states_list))
    if ((is.null(inputs$list_of_distances_mats) == FALSE)) {
      distances_mat = inputs$list_of_distances_mats[[i]]
    }
    else {
      distances_mat = matrix(1, nrow = length(areas), 
                             ncol = length(areas))
    }
    x = BioGeoBEARS_model_object@params_table["x", "est"]
    dispersal_multipliers_matrix = distances_mat^x
    if ((is.null(inputs$list_of_envdistances_mats) == FALSE)) {
      envdistances_mat = inputs$list_of_envdistances_mats[[i]]
    }
    else {
      envdistances_mat = matrix(1, nrow = length(areas), 
                                ncol = length(areas))
    }
    n = BioGeoBEARS_model_object@params_table["n", "est"]
    dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
      envdistances_mat^n
    if ((is.null(inputs$list_of_dispersal_multipliers_mats) == 
         FALSE)) {
      manual_dispersal_multipliers_matrix = as.matrix(inputs$list_of_dispersal_multipliers_mats[[i]])
    }
    else {
      manual_dispersal_multipliers_matrix = matrix(1, 
                                                   nrow = length(areas), ncol = length(areas))
    }
    w = BioGeoBEARS_model_object@params_table["w", "est"]
    dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
      manual_dispersal_multipliers_matrix^w
    dmat_times_d = dispersal_multipliers_matrix * matrix(d, 
                                                         nrow = length(areas), ncol = length(areas))
    amat = dispersal_multipliers_matrix * matrix(a, nrow = length(areas), 
                                                 ncol = length(areas))
    if ((is.null(inputs$list_of_area_of_areas) == FALSE)) {
      area_of_areas = inputs$list_of_area_of_areas[[i]]
    }
    else {
      area_of_areas = rep(1, length(areas))
    }
    u = BioGeoBEARS_model_object@params_table["u", "est"]
    extinction_modifier_list = area_of_areas^(1 * u)
    elist = extinction_modifier_list * rep(e, length(areas))
    if (traitTF == FALSE) {
      Qmat_tmp = rcpp_states_list_to_DEmat(areas_list = allareas_list, 
                                           states_list = all_states_list[states_allowed_TF], 
                                           dmat = dmat_times_d, elist = elist, amat = amat, 
                                           include_null_range = include_null_range, normalize_TF = TRUE, 
                                           makeCOO_TF = force_sparse)
    }
    if (traitTF == TRUE) {
      num_geog_states = length(all_states_list[states_allowed_TF])
      numstates_geogtrait = num_trait_states * num_geog_states
      tmpres = modify_Qmat_with_trait(Qmat = NULL, BioGeoBEARS_run_object, 
                                      numstates_geogtrait = numstates_geogtrait, areas_list = allareas_list, 
                                      states_list = all_states_list[states_allowed_TF], 
                                      dispersal_multipliers_matrix = dispersal_multipliers_matrix, 
                                      elist = elist, force_sparse = force_sparse)
      Qmat_tmp = tmpres$Qmat
      m = tmpres$m
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
    else {
      num_geog_states = length(all_states_list[states_allowed_TF])
      numstates_geogtrait = num_geog_states
    }
    if (force_sparse == TRUE) {
      tmpQmat_in_REXPOKIT_coo_fmt = Qmat_tmp
      tmpQmat_in_kexpmv_crs_fmt = coo2crs(ia = tmpQmat_in_REXPOKIT_coo_fmt[, 
                                                                           "ia"], ja = tmpQmat_in_REXPOKIT_coo_fmt[, "ja"], 
                                          a = tmpQmat_in_REXPOKIT_coo_fmt[, "a"], n = numstates_geogtrait, 
                                          transpose_needed = FALSE)
      Qmat_tmp = tmpQmat_in_REXPOKIT_coo_fmt
    }
    if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 
        1) {
      tr = check_trfn(trfn = inputs$trfn)
      tree_to_chainsaw = NULL
      tree_to_chainsaw[[1]] = tr
      return_pieces_list = NULL
      return_pieces_list[[1]] = tr
      return_pieces_basenames = NULL
      tmp_labels_merge = paste(tr$tip.label, collapse = ",", 
                               sep = "")
      tmp_labels_split = strsplit(tmp_labels_merge, split = ",")[[1]]
      return_pieces_basenames[[1]] = paste(sort(tmp_labels_split), 
                                           collapse = ",", sep = "")
      chainsaw_object = list()
      chainsaw_object$tree_to_chainsaw = tree_to_chainsaw
      chainsaw_object$return_pieces_list = return_pieces_list
      chainsaw_object$return_pieces_basenames = return_pieces_basenames
      attr(chainsaw_object, "class") = "chainsaw_result"
      inputs$tree_sections_list[[1]] = chainsaw_object
    }
    spPmat_inputs = get_spPmat_inputs_from_BGB(BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                               states_list = all_states_list[states_allowed_TF], 
                                               dispersal_multipliers_matrix = dispersal_multipliers_matrix)
    chainsaw_result = inputs$tree_sections_list[[i]]
    current_tip_relative_probs_of_each_state
    new_tip_likelihoods = matrix(0, nrow = length(chainsaw_result$return_pieces_list), 
                                 ncol = ncol(current_tip_relative_probs_of_each_state))
    if (traitTF == TRUE) {
      wTrait_states_allowed_TF = c(rep(states_allowed_TF, 
                                       times = num_trait_states))
      if (sum(wTrait_states_allowed_TF) != numstates_geogtrait) {
        txt = paste0("STOP ERROR in calc_loglike_sp_stratified(): sum(wTrait_states_allowed_TF)=", 
                     sum(wTrait_states_allowed_TF), ", and numstates_geogtrait=", 
                     numstates_geogtrait, ". They must be equal to proceed.")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
    }
    for (jj in 1:length(chainsaw_result$return_pieces_list)) {
      treepiece = chainsaw_result$return_pieces_list[[jj]]
      treepiece
      if (is.numeric(treepiece)) {
        do_exponentiation = TRUE
        if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == 
                                                 TRUE)) {
          TF1 = inputs$master_table$stratum == i
          TF2 = inputs$master_table$piecenum == jj
          TF3 = inputs$master_table$piececlass == "subbranch"
          TF = (TF1 + TF2 + TF3) == 3
          this_row_of_master_table_is_being_used = TF
          rownum = (1:nrow(condlikes_table))[TF]
          tmp_master_table_row = inputs$master_table[rownum, 
          ]
          if (nrow(tmp_master_table_row) != 1) {
            stoptxt = paste("\n\nFATAL ERROR in stratified loglike calculation at i=", 
                            i, "; jj=", jj, "; ", "inputs$master_table$piececlass == \"subbranch\"", 
                            "\nnrow(tmp_master_table_row) should =1 but instead =", 
                            nrow(tmp_master_table_row), "\n", sep = "")
            stop(stoptxt)
          }
          master_tip_time_bp = tmp_master_table_row$time_bp
          time_top = tmp_master_table_row$time_top
          time_bot = tmp_master_table_row$time_bot
          is_fossil = tmp_master_table_row$fossils
          if ((master_tip_time_bp > time_top) && (is.na(is_fossil) == 
                                                  FALSE) && (is_fossil == TRUE)) {
            do_exponentiation = FALSE
          }
          if ((master_tip_time_bp >= time_top) && (master_tip_time_bp < 
                                                   time_bot) && (is.na(is_fossil) == FALSE) && 
              (is_fossil == TRUE)) {
            amount_to_shorten_by = master_tip_time_bp - 
              time_top
            treepiece = treepiece - amount_to_shorten_by
            do_exponentiation = TRUE
          }
          if (tmp_master_table_row$edge.length < min_branchlength) {
            do_exponentiation = FALSE
          }
        }
        tipname = chainsaw_result$return_pieces_basenames[[jj]]
        tip_TF = phy_as_it_is_chopped_down$tip.label == 
          tipname
        relative_probs_of_each_state_at_the_tip_of_this_branch = current_tip_relative_probs_of_each_state[tip_TF, 
                                                                                                          states_to_use_TF]
        if (do_exponentiation == TRUE) {
          if (force_sparse == FALSE) {
            independent_likelihoods_at_branch_section_bottom = expokit_dgpadm_Qmat2(times = treepiece, 
                                                                                    Qmat = Qmat_tmp, transpose_needed = TRUE)
            if (traitTF == FALSE) {
              tmp_conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% 
                                                                              relative_probs_of_each_state_at_the_tip_of_this_branch[states_allowed_TF], 
                                                                            nrow = 1)
            }
            else {
              tmp_conditional_likelihoods_at_branch_section_bottom = matrix(independent_likelihoods_at_branch_section_bottom %*% 
                                                                              relative_probs_of_each_state_at_the_tip_of_this_branch[wTrait_states_allowed_TF], 
                                                                            nrow = 1)
            }
          }
          else {
            if (class(Qmat_tmp) != "data.frame") {
              txt = paste0("ERROR: calc_loglike_sp_stratified is attempting to use a sparse COO-formated Q matrix, to calculated the likelihoods down one branch segment, but you provided a Qmat not in data.frame form")
              cat("\n\n")
              cat(txt)
              cat("\n\n")
              print("class(Qmat_tmp)")
              print(class(Qmat_tmp))
              print("Qmat_tmp")
              print(Qmat_tmp)
              stop(txt)
            }
            if ((ncol(Qmat_tmp) != 3)) {
              txt = paste0("ERROR: calc_loglike_sp_stratified is attempting to use a sparse COO-formated Q matrix, to calculated the likelihoods down one branch segment, but you provided a Qmat that does't have 3 columns")
              cat("\n\n")
              cat(txt)
              cat("\n\n")
              print("class(Qmat_tmp)")
              print(class(Qmat_tmp))
              print("Qmat_tmp")
              print(Qmat_tmp)
              stop(txt)
            }
            coo_n = numstates_geogtrait
            anorm = 1
            if (traitTF == FALSE) {
              try_result_segment = try(kexpmv::expokit_dgexpv(mat = tmpQmat_in_kexpmv_crs_fmt, 
                                                              t = treepiece, vector = relative_probs_of_each_state_at_the_tip_of_this_branch[states_allowed_TF], 
                                                              transpose_needed = TRUE, transform_to_crs = FALSE, 
                                                              crs_n = numstates_geogtrait, anorm = NULL, 
                                                              mxstep = 10000, tol = as.numeric(1e-10)))
            }
            else {
              try_result_segment = try(kexpmv::expokit_dgexpv(mat = tmpQmat_in_kexpmv_crs_fmt, 
                                                              t = treepiece, vector = relative_probs_of_each_state_at_the_tip_of_this_branch[wTrait_states_allowed_TF], 
                                                              transpose_needed = TRUE, transform_to_crs = FALSE, 
                                                              crs_n = numstates_geogtrait, anorm = NULL, 
                                                              mxstep = 10000, tol = as.numeric(1e-10)))
            }
            if (printlevel >= 1) {
              txt = "S"
              cat(txt)
            }
            if (class(try_result_segment) == "try-error") {
              cat("\n\ntry-error on kexpmv::expokit_dgexpv():\n\n")
              cat("i=", i, "\n")
              cat("t=treepiece==", treepiece, "\n")
              print(tmpQmat_in_kexpmv_crs_fmt)
              print(phy2)
              print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
              print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum, 
              ])
              print(coo_n)
              print(anorm)
              print("BioGeoBEARS_model_object")
              print(BioGeoBEARS_model_object)
              stoptxt = "\n\nStopping on error in sparse exponentiation downpass (treepiece, aka a branch segment): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nAnother solution: To have this error report an extremely low log-likelihood,, set BioGeoBEARS_run_object$on_NaN_error to something like -1e50.\n\n"
              if (is.null(on_NaN_error)) {
                stop(stoptxt)
              }
              print("print(on_NaN_error):")
              print(on_NaN_error)
              if ((is.numeric(on_NaN_error)) && (return_what == 
                                                 "loglike")) {
                warning(paste0("\n\nWarning  on error in sparse exponentiation downpass (treepiece, aka a branch segment): NaNs produced in likelihood calculation. This may mean your transition matrix disallows necessary transitions.  E.g., if your ranges are 'A' and 'B', and your model is DEC, then allowing range 'AB' as a possible state is required, so that you can get from 'A' to 'B' via 'AB' as the intermediate. Alternatively, NaNs can be produced sometimes if your Maximum Likelihood (ML) search proposes weird parameter values (such as a negative rate or weight) or a parameter so small that required transitions have a probability that machine precision rounds to zero or negative.  Sometimes this seems to occur because optim, optimx, etc. propose parameters slightly outside the user-specified upper and lower (min/max) boundaries for some reason. One solution is often to narrow the min/max limits. \n\nYou are using another solution: Normally, this would be a stop error, but you specified that BioGeoBEARS_run_object$on_NaN_error=", 
                               on_NaN_error, "\n\n"))
                return(on_NaN_error)
              }
              else {
                stop(stoptxt)
              }
            }
            tmp_conditional_likelihoods_at_branch_section_bottom = c(try_result_segment$output_probs)
            tmp_conditional_likelihoods_at_branch_section_bottom[tmp_conditional_likelihoods_at_branch_section_bottom < 
                                                                   0] = tmp_conditional_likelihoods_at_branch_section_bottom
          }
          if (traitTF == FALSE) {
            conditional_likelihoods_at_branch_section_bottom = matrix(0, 
                                                                      nrow = 1, ncol = length(relative_probs_of_each_state_at_the_tip_of_this_branch))
            conditional_likelihoods_at_branch_section_bottom[states_allowed_TF] = tmp_conditional_likelihoods_at_branch_section_bottom
          }
          else {
            conditional_likelihoods_at_branch_section_bottom = matrix(0, 
                                                                      nrow = 1, ncol = length(relative_probs_of_each_state_at_the_tip_of_this_branch))
            conditional_likelihoods_at_branch_section_bottom[wTrait_states_allowed_TF] = tmp_conditional_likelihoods_at_branch_section_bottom
          }
          conditional_likelihoods_at_branch_section_bottom[states_allowed_TF == 
                                                             FALSE] = 0
        }
        else {
          conditional_likelihoods_at_branch_section_bottom = matrix(relative_probs_of_each_state_at_the_tip_of_this_branch, 
                                                                    nrow = 1)
        }
        chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom
        chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]] = conditional_likelihoods_at_branch_section_bottom/sum(conditional_likelihoods_at_branch_section_bottom)
        chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]] = conditional_likelihoods_at_branch_section_bottom/sum(conditional_likelihoods_at_branch_section_bottom)
        if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == 
                                                 TRUE)) {
          TF1 = inputs$master_table$stratum == i
          TF2 = inputs$master_table$piecenum == jj
          TF3 = inputs$master_table$piececlass == "subbranch"
          TF4 = inputs$master_table$piececlass == "orig_tip"
          TF5 = (TF3 + TF4) == 1
          TF = (TF1 + TF2 + TF5) == 3
          rownum = (1:nrow(condlikes_table))[TF]
          condlikes_table[rownum, ] = conditional_likelihoods_at_branch_section_bottom
          if (calc_ancprobs == TRUE) {
            relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[rownum, 
            ] = conditional_likelihoods_at_branch_section_bottom/sum(conditional_likelihoods_at_branch_section_bottom)
          }
        }
      }
      else {
        tmp_subtree = treepiece
        if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == 
                                                 TRUE)) {
          tmp_subtree_tipnums = 1:length(tmp_subtree$tip.label)
          for (iter in 1:length(tmp_subtree_tipnums)) {
            subtree_tip = tmp_subtree_tipnums[iter]
            TF1 = inputs$master_table$stratum == i
            TF2 = inputs$master_table$piecenum == jj
            TF3 = inputs$master_table$piececlass == 
              "subtree"
            TF4 = inputs$master_table$SUBnode == subtree_tip
            TF = (TF1 + TF2 + TF3 + TF4) == 4
            this_row_of_master_table_is_being_used = TF
            rownum = (1:nrow(inputs$master_table))[TF]
            tmp_master_table_row = inputs$master_table[rownum, 
            ]
            if (nrow(tmp_master_table_row) != 1) {
              stoptxt = paste("\n\nFATAL ERROR in stratified loglike calculation at i=", 
                              i, "; jj=", jj, "; ", "inputs$master_table$piececlass == \"subtree\"", 
                              "; subtree_tip=", subtree_tip, "\nnrow(tmp_master_table_row) should =1 but instead =", 
                              nrow(tmp_master_table_row), "\n", sep = "")
              stop(stoptxt)
            }
            master_tip_time_bp = tmp_master_table_row$time_bp
            time_top = tmp_master_table_row$time_top
            time_bot = tmp_master_table_row$time_bot
            is_fossil = tmp_master_table_row$fossils
            if ((master_tip_time_bp >= time_top) && 
                (master_tip_time_bp < time_bot) && (is_fossil == 
                                                    TRUE)) {
              amount_to_shorten_by = master_tip_time_bp - 
                time_top
              tmp2_edgeTF = tmp_subtree$edge[, 2] == 
                subtree_tip
              tmp2_edgenum = (1:nrow(tmp_subtree$edge))[tmp2_edgeTF]
            }
          }
        }
        tipnames = tmp_subtree$tip.label
        tips_for_subtree_TF = phy_as_it_is_chopped_down$tip.label %in% 
          tipnames
        if (traitTF == FALSE) {
          subtree_tip_relative_probs_of_each_state = current_tip_relative_probs_of_each_state[tips_for_subtree_TF, 
          ][, states_allowed_TF]
        }
        else {
          subtree_tip_relative_probs_of_each_state = current_tip_relative_probs_of_each_state[tips_for_subtree_TF, 
          ][, wTrait_states_allowed_TF]
        }
        if (sum(states_allowed_TF) == 1) {
          if (traitTF == FALSE) {
            subtree_tip_relative_probs_of_each_state = matrix(data = subtree_tip_relative_probs_of_each_state, 
                                                              ncol = 1)
          }
          else {
            subtree_tip_relative_probs_of_each_state = matrix(data = subtree_tip_relative_probs_of_each_state, 
                                                              ncol = sum(wTrait_states_allowed_TF))
          }
        }
        tmp_fixnode = NULL
        tmp_fixlikes = NULL
        if ((!is.null(fixnode)) && (length(fixnode) > 
                                    0)) {
          if (length(fixnode) > 1) {
            TF1 = inputs$master_table$stratum == i
            TF2 = inputs$master_table$piecenum == jj
            TF3 = inputs$master_table$piececlass == 
              "subtree"
            TF = ((TF1 + TF2 + TF3) == 3)
            tmprows = inputs$master_table[TF, ]
            fixnodes_in_subtree_TF = fixnode %in% tmprows$node
            if (sum(fixnodes_in_subtree_TF) > 0) {
              temporary_fixnodes = fixnode[fixnodes_in_subtree_TF]
              if (traitTF == FALSE) {
                temporary_fixlikes = fixlikes[fixnodes_in_subtree_TF, 
                                              states_allowed_TF]
              }
              else {
                temporary_fixlikes = fixlikes[fixnodes_in_subtree_TF, 
                                              wTrait_states_allowed_TF]
              }
              subtree_rows_in_fixnodes_TF = tmprows$node %in% 
                fixnode
              subtree_fixnode_master_nodenums = tmprows$node[subtree_rows_in_fixnodes_TF]
              subtree_fixnode_nums = tmprows$SUBnode[subtree_rows_in_fixnodes_TF]
              order_subtree_fixnode_nums = order(subtree_fixnode_nums)
              subtree_fixnode_nums = subtree_fixnode_nums[order_subtree_fixnode_nums]
              if (length(order_subtree_fixnode_nums) > 
                  1) {
                temporary_fixlikes = temporary_fixlikes[order_subtree_fixnode_nums, 
                ]
              }
            }
            else {
              temporary_fixnodes = NULL
              subtree_fixnode_master_nodenums = NULL
              subtree_fixnode_nums = NULL
              temporary_fixlikes = NULL
            }
            TF1 = unique(tmprows$stratum) == i
            TF2 = unique(tmprows$piecenum) == jj
            TF3 = unique(tmprows$piececlass) == "subtree"
            TF = ((TF1 + TF2 + TF3) == 3)
            if (TF == TRUE) {
              if (length(subtree_fixnode_nums) == 0) {
                subtree_fixnode_nums = NULL
                temporary_fixlikes = NULL
              }
              tmp_fixnode = subtree_fixnode_nums
              tmp_fixlikes = temporary_fixlikes
            }
            else {
              tmp_fixnode = NULL
              tmp_fixlikes = NULL
            }
          }
          else {
            temporary_fixnode = fixnode
            temporary_fixlikes = c(fixlikes)
            TF1 = inputs$master_table$node == temporary_fixnode
            TF2 = inputs$master_table$SUBnode.type == 
              "root"
            TF3 = inputs$master_table$SUBnode.type == 
              "internal"
            TF = ((TF1 + TF2 + TF3) == 2)
            tmprow = inputs$master_table[TF, ]
            TF1 = tmprow$stratum == i
            TF2 = tmprow$piecenum == jj
            TF3 = tmprow$piececlass == "subtree"
            TF = ((TF1 + TF2 + TF3) == 3)
            if (TF == TRUE) {
              tmp_fixnode = tmprow$SUBnode
              tmp_fixlikes = temporary_fixlikes[states_allowed_TF]
            }
            else {
              tmp_fixnode = NULL
              tmp_fixlikes = NULL
            }
          }
        }
        calc_loglike_sp_results = calc_loglike_sp(tip_condlikes_of_data_on_each_state = subtree_tip_relative_probs_of_each_state, 
                                                  phy = tmp_subtree, Qmat = Qmat_tmp, spPmat = NULL, 
                                                  min_branchlength = min_branchlength, return_what = "all", 
                                                  probs_of_states_at_root = NULL, rootedge = TRUE, 
                                                  sparse = force_sparse, printlevel = printlevel, 
                                                  use_cpp = TRUE, input_is_COO = force_sparse, 
                                                  spPmat_inputs = spPmat_inputs, cppSpMethod = cppSpMethod, 
                                                  cluster_already_open = cluster_already_open, 
                                                  calc_ancprobs = calc_ancprobs, include_null_range = include_null_range, 
                                                  fixnode = tmp_fixnode, fixlikes = tmp_fixlikes, 
                                                  stratified = TRUE, states_allowed_TF = rep(TRUE, 
                                                                                             times = ncol(subtree_tip_relative_probs_of_each_state)), 
                                                  m = m, jts_matrix = jts_matrix, BioGeoBEARS_model_object = BioGeoBEARS_model_object, 
                                                  on_NaN_error = BioGeoBEARS_run_object$on_NaN_error)
        if (!is.null(inputs$lists_of_states_lists_0based)) {
          names_of_calc_loglike_sp_results_objects = names(calc_loglike_sp_results)
          for (name_i in 1:length(calc_loglike_sp_results)) {
            oldmat = calc_loglike_sp_results[[name_i]]
            TF1 = names_of_calc_loglike_sp_results_objects[name_i] == 
              "relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS"
            TF2 = names_of_calc_loglike_sp_results_objects[name_i] == 
              "condlikes_of_each_state"
            TF3 = names_of_calc_loglike_sp_results_objects[name_i] == 
              "relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS"
            TF4 = names_of_calc_loglike_sp_results_objects[name_i] == 
              "relative_probs_of_each_state_at_bottom_of_root_branch"
            if (TF1 || TF2 || TF3) {
              if (traitTF == FALSE) {
                newmat = matrix(0, nrow = nrow(oldmat), 
                                ncol = length(states_allowed_TF))
                newmat[, states_allowed_TF] = oldmat
              }
              if (traitTF == TRUE) {
                full_matrix_ncols = length(states_allowed_TF) * 
                  num_trait_states
                newmat = matrix(0, nrow = nrow(oldmat), 
                                ncol = full_matrix_ncols)
                wTrait_states_allowed_TF = c(rep(states_allowed_TF, 
                                                 times = num_trait_states))
                newmat[, wTrait_states_allowed_TF] = oldmat
              }
              calc_loglike_sp_results[[name_i]] = newmat
            }
            if (TF4) {
              if (traitTF == FALSE) {
                newmat = matrix(0, nrow = 1, ncol = length(states_allowed_TF))
                newmat[, states_allowed_TF] = oldmat
              }
              if (traitTF == TRUE) {
                full_matrix_ncols = length(states_allowed_TF) * 
                  num_trait_states
                newmat = matrix(0, nrow = 1, ncol = full_matrix_ncols)
                wTrait_states_allowed_TF = c(rep(states_allowed_TF, 
                                                 times = num_trait_states))
                newmat[, wTrait_states_allowed_TF] = oldmat
              }
              calc_loglike_sp_results[[name_i]] = newmat
            }
          }
        }
        tmp_tipnums = 1:length(tipnames)
        if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == 
                                                 TRUE)) {
          for (rownum in 1:nrow(calc_loglike_sp_results$condlikes_of_each_state)) {
            tmp_condlikes = calc_loglike_sp_results$condlikes_of_each_state[rownum, 
            ]
            subtree_node = rownum
            TF1 = inputs$master_table$stratum == i
            TF2 = inputs$master_table$piecenum == jj
            TF3 = inputs$master_table$piececlass == 
              "subtree"
            TF4 = inputs$master_table$SUBnode == subtree_node
            TF = (TF1 + TF2 + TF3 + TF4) == 4
            condlikes_table_rownum = (1:nrow(condlikes_table))[TF]
            condlikes_table[condlikes_table_rownum, 
            ] = tmp_condlikes
            if (calc_ancprobs == TRUE) {
              if (rownum <= nrow(calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)) {
                if (inputs$master_table$SUBnode.type[condlikes_table_rownum] != 
                    "root") {
                  tmp = calc_loglike_sp_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[rownum, 
                  ]
                  relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[condlikes_table_rownum, 
                  ] = tmp
                }
                else {
                  tmp = calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch
                  relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[condlikes_table_rownum, 
                  ] = tmp
                }
              }
            }
          }
        }
        chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]] = matrix(data = calc_loglike_sp_results$condlikes_of_each_state[-tmp_tipnums, 
        ], ncol = ncol(calc_loglike_sp_results$condlikes_of_each_state))
        chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]] = calc_loglike_sp_results$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[-tmp_tipnums, 
        ]
        chainsaw_result$relative_probs_of_each_state_at_bottom_of_root_branch[[jj]] = calc_loglike_sp_results$relative_probs_of_each_state_at_bottom_of_root_branch
        numrows_to_add = nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
        rownum_for_bottom_of_root = nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]])
        startrow = current_condlikes_row + 1
        endrow = current_condlikes_row + numrows_to_add
        all_relative_probs_of_each_state[startrow:endrow, 
                                         states_to_use_TF] = chainsaw_result$relative_probabilities_for_nodes_plus_bottom_in_this_section[[jj]]
        all_condlikes_of_each_state[startrow:endrow, 
                                    states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]
        current_condlikes_row = current_condlikes_row + 
          numrows_to_add
      }
      new_tip_likelihoods[jj, states_to_use_TF] = chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]][nrow(chainsaw_result$conditional_likelihoods_for_nodes_plus_bottom_in_this_section[[jj]]), 
      ]
    }
    current_tip_relative_probs_of_each_state = new_tip_likelihoods
    phy_as_it_is_chopped_down = chainsaw_result$tree_to_chainsaw
  }
  all_condlikes_of_each_state_zero_TF = all_condlikes_of_each_state == 
    0
  all_condlikes_of_each_state_nonzero_TF = all_condlikes_of_each_state_zero_TF == 
    FALSE
  rows_that_are_NOT_numeric_zeros_TF = rowSums(all_condlikes_of_each_state_nonzero_TF) >= 
    1
  final_all_condlikes_of_each_state = all_condlikes_of_each_state[rows_that_are_NOT_numeric_zeros_TF, 
  ]
  all_relative_probs_of_each_state = all_relative_probs_of_each_state[rows_that_are_NOT_numeric_zeros_TF, 
  ]
  if (rootedge == TRUE) {
    grand_total_likelihood = sum(log(rowSums(final_all_condlikes_of_each_state)))
    grand_total_likelihood
  }
  else {
    grand_total_likelihood = sum(log(rowSums(final_all_condlikes_of_each_state[-nrow(final_all_condlikes_of_each_state), 
    ])))
    grand_total_likelihood
  }
  if (is.na(grand_total_likelihood) == TRUE) {
    TF = is.na(all_relative_probs_of_each_state[, 1])
    tmpr = (1:nrow(all_relative_probs_of_each_state))[TF]
    stoptxt1 = paste("\n\nFATAL ERROR IN calc_loglike_sp_stratified(). grand_total_likelihood=NA.\n", 
                     "These rows of 'all_relative_probs_of_each_state' had NAs:\n", 
                     paste(tmpr, collpase = ",", sep = ""), "\n", "\n", 
                     "One possible cause of this: your dispersal matrix may be too restrictive; try changing\n", 
                     "e.g. the 0 values to e.g. 0.0000001.  Good luck!", 
                     sep = "")
    if (printlevel > 0) {
      cat(stoptxt1)
    }
    if (is.null(BioGeoBEARS_run_object$on_NaN_error) == 
        TRUE) {
      stop(stoptxt1)
    }
    else {
      grand_total_likelihood = BioGeoBEARS_run_object$on_NaN_error
    }
  }
  if (calc_TTL_loglike_from_condlikes_table == TRUE) {
    TF2 = inputs$master_table$SUBnode.type == "internal"
    TF3 = inputs$master_table$SUBnode.type == "orig_tip"
    TF4 = inputs$master_table$SUBnode.type == "root"
    TF234 = (TF2 + TF3 + TF4) == 1
    sum(TF234)
    TF = TF234 == 1
    nodes_in_original_tree = inputs$master_table[TF, ]
    node_order_original = order(nodes_in_original_tree$node)
    condlikes_of_each_state = condlikes_table[TF, ][node_order_original, 
    ]
    computed_likelihoods_at_each_node = rowSums(condlikes_of_each_state)
    grand_total_likelihood = sum(log(computed_likelihoods_at_each_node))
    if (calc_ancprobs == TRUE) {
      TF2 = inputs$master_table$SUBnode.type == "internal"
      TF3 = inputs$master_table$SUBnode.type == "tip"
      TF4 = inputs$master_table$piececlass == "subtree"
      TF234 = (TF2 + TF3 + TF4) == 2
      sum(TF234)
      TF = TF234 == 1
      nodes_in_original_tree = inputs$master_table[TF, 
      ]
      node_order_original = order(nodes_in_original_tree$node)
      tmptable = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE[TF, 
      ][node_order_original, ]
      relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = tmptable/rowSums(tmptable)
      root_row = rep(NA, times = ncol(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
      tmpmat1 = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[1:length(original_phy$tip.label), 
      ]
      tmpmat3_rows = (length(original_phy$tip.label) + 
                        1):nrow(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)
      tmpmat3 = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[tmpmat3_rows, 
      ]
      relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = rbind(tmpmat1, 
                                                                                root_row, tmpmat3)
      relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = NULL
      tmptable = condlikes_of_each_state/rowSums(condlikes_of_each_state)
      relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = tmptable
      anc_row_of_master_table_TF = inputs$master_table$node.type == 
        "root"
      anc_node_original_tree = inputs$master_table$node[anc_row_of_master_table_TF]
      anc_node_original_tree
      starting_probs = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc_node_original_tree, 
      ]
    }
  }
  if (calc_ancprobs == TRUE) {
    cat("\nUppass started for (STRATIFIED) marginal ancestral states estimation!\n", 
        sep = "")
    numrows_for_UPPASS = original_phy$Nnode + length(original_phy$tip.label)
    relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = matrix(data = 0, 
                                                                       nrow = numrows_for_UPPASS, ncol = ncol(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS))
    relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = matrix(data = 0, 
                                                                             nrow = numrows_for_UPPASS, ncol = ncol(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS))
    starting_probs
    relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, 
    ] = 1/length(starting_probs)
    for (i in num_timeperiods:1) {
      user_specified_constraints_on_states_list_TF = FALSE
      states_allowed_TF1 = rep(TRUE, length(all_states_list))
      states_allowed_TF2 = rep(TRUE, length(all_states_list))
      states_allowed_TF3 = rep(TRUE, length(all_states_list))
      if ((is.null(inputs$list_of_areas_allowed_mats) == 
           FALSE)) {
        user_specified_constraints_on_states_list_TF = TRUE
      }
      if ((is.null(inputs$list_of_areas_adjacency_mats) == 
           FALSE)) {
        user_specified_constraints_on_states_list_TF = TRUE
      }
      if ((is.null(inputs$lists_of_states_lists_0based) == 
           FALSE)) {
        user_specified_constraints_on_states_list_TF = TRUE
      }
      if (user_specified_constraints_on_states_list_TF == 
          TRUE) {
        if ((is.null(inputs$list_of_areas_allowed_mats) == 
             FALSE)) {
          areas_allowed_mat = inputs$list_of_areas_allowed_mats[[i]]
          cat("\ni=", i, "\n", sep = "")
          cat("areas_allowed_mat: ", sep = "")
          print(areas_allowed_mat)
          states_allowed_TF1 = sapply(X = all_states_list, 
                                      FUN = check_if_state_is_allowed, areas_allowed_mat)
          if (include_null_range == TRUE) {
            states_allowed_TF1[1] = TRUE
          }
        }
        if ((is.null(inputs$list_of_areas_adjacency_mats) == 
             FALSE)) {
          areas_adjacency_mat = inputs$list_of_areas_adjacency_mats[[i]]
          states_allowed_TF2 = sapply(X = all_states_list, 
                                      FUN = check_if_state_is_allowed_by_adjacency, 
                                      areas_adjacency_mat)
          if (include_null_range == TRUE) {
            states_allowed_TF2[1] = TRUE
          }
        }
        if ((is.null(inputs$lists_of_states_lists_0based) == 
             FALSE)) {
          states_allowed_TF3 = all_states_list %in% 
            inputs$lists_of_states_lists_0based[[i]]
          if (include_null_range == TRUE) {
            states_allowed_TF3[1] = TRUE
          }
        }
        states_allowed_TF = ((states_allowed_TF1 + states_allowed_TF2 + 
                                states_allowed_TF3) == 3)
        inputs$lists_of_states_lists_0based[[i]] = all_states_list[states_allowed_TF]
      }
      else {
        pass = 1
        states_allowed_TF = rep(TRUE, length(all_states_list))
      }
      states_to_use_TF = rep(TRUE, length(all_states_list))
      if ((is.null(inputs$list_of_distances_mats) == FALSE)) {
        distances_mat = inputs$list_of_distances_mats[[i]]
      }
      else {
        distances_mat = matrix(1, nrow = length(areas), 
                               ncol = length(areas))
      }
      dispersal_multipliers_matrix = distances_mat^x
      if ((is.null(inputs$list_of_envdistances_mats) == 
           FALSE)) {
        envdistances_mat = inputs$list_of_envdistances_mats[[1]]
      }
      else {
        envdistances_mat = matrix(1, nrow = length(areas), 
                                  ncol = length(areas))
      }
      n = BioGeoBEARS_model_object@params_table["n", "est"]
      dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
        envdistances_mat^n
      if ((is.null(inputs$list_of_dispersal_multipliers_mats) == 
           FALSE)) {
        manual_dispersal_multipliers_matrix = as.matrix(inputs$list_of_dispersal_multipliers_mats[[i]])
      }
      else {
        manual_dispersal_multipliers_matrix = matrix(1, 
                                                     nrow = length(areas), ncol = length(areas))
      }
      w = BioGeoBEARS_model_object@params_table["w", "est"]
      dispersal_multipliers_matrix = dispersal_multipliers_matrix * 
        manual_dispersal_multipliers_matrix^w
      dmat_times_d = dispersal_multipliers_matrix * matrix(d, 
                                                           nrow = length(areas), ncol = length(areas))
      amat = dispersal_multipliers_matrix * matrix(a, 
                                                   nrow = length(areas), ncol = length(areas))
      if ((is.null(inputs$list_of_area_of_areas) == FALSE)) {
        area_of_areas = inputs$list_of_area_of_areas[[i]]
      }
      else {
        area_of_areas = rep(1, length(areas))
      }
      extinction_modifier_list = area_of_areas^(1 * u)
      elist = extinction_modifier_list * rep(e, length(areas))
      if (traitTF == FALSE) {
        Qmat_tmp = rcpp_states_list_to_DEmat(areas_list = allareas_list, 
                                             states_list = all_states_list[states_allowed_TF], 
                                             dmat = dmat_times_d, elist = elist, amat = amat, 
                                             include_null_range = include_null_range, normalize_TF = TRUE, 
                                             makeCOO_TF = force_sparse)
      }
      if (traitTF == TRUE) {
        num_geog_states = length(all_states_list[states_allowed_TF])
        numstates_geogtrait = num_trait_states * num_geog_states
        wTrait_states_allowed_TF = c(rep(states_allowed_TF, 
                                         times = num_trait_states))
        if (ncol(tip_condlikes_of_data_on_each_state[, 
                                                     wTrait_states_allowed_TF]) != numstates_geogtrait) {
          txt = paste0("STOP ERROR in calc_loglike_sp_stratified(): ncol(tip_condlikes_of_data_on_each_state)=", 
                       ncol(tip_condlikes_of_data_on_each_state), 
                       ", and numstates_geogtrait=", numstates_geogtrait, 
                       ". They must be equal to proceed.")
          cat("\n\n")
          cat(txt)
          cat("\n\n")
          stop(txt)
        }
        tmpres = modify_Qmat_with_trait(Qmat = NULL, 
                                        BioGeoBEARS_run_object, numstates_geogtrait = numstates_geogtrait, 
                                        areas_list = allareas_list, states_list = all_states_list[states_allowed_TF], 
                                        dispersal_multipliers_matrix = dispersal_multipliers_matrix, 
                                        elist = elist, force_sparse = force_sparse)
        Qmat_tmp = tmpres$Qmat
        m = tmpres$m
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
      else {
        num_geog_states = length(all_states_list[states_allowed_TF])
        numstates_geogtrait = num_geog_states
      }
      if (force_sparse == TRUE) {
        tmpQmat_in_REXPOKIT_coo_fmt = Qmat_tmp
        tmpQmat_in_kexpmv_crs_fmt = coo2crs(ia = tmpQmat_in_REXPOKIT_coo_fmt[, 
                                                                             "ia"], ja = tmpQmat_in_REXPOKIT_coo_fmt[, 
                                                                                                                     "ja"], a = tmpQmat_in_REXPOKIT_coo_fmt[, "a"], 
                                            n = numstates_geogtrait, transpose_needed = FALSE)
      }
      if (is.null(inputs$timeperiods) || length(inputs$timeperiods) == 
          1) {
        tr = check_trfn(trfn = inputs$trfn)
        tree_to_chainsaw = NULL
        tree_to_chainsaw[[1]] = tr
        return_pieces_list = NULL
        return_pieces_list[[1]] = tr
        return_pieces_basenames = NULL
        tmp_labels_merge = paste(tr$tip.label, collapse = ",", 
                                 sep = "")
        tmp_labels_split = strsplit(tmp_labels_merge, 
                                    split = ",")[[1]]
        return_pieces_basenames[[1]] = paste(sort(tmp_labels_split), 
                                             collapse = ",", sep = "")
        chainsaw_object = list()
        chainsaw_object$tree_to_chainsaw = tree_to_chainsaw
        chainsaw_object$return_pieces_list = return_pieces_list
        chainsaw_object$return_pieces_basenames = return_pieces_basenames
        attr(chainsaw_object, "class") = "chainsaw_result"
        inputs$tree_sections_list[[1]] = chainsaw_object
      }
      spPmat_inputs = get_spPmat_inputs_from_BGB(BioGeoBEARS_run_object = BioGeoBEARS_run_object, 
                                                 states_list = all_states_list[states_allowed_TF], 
                                                 dispersal_multipliers_matrix = dispersal_multipliers_matrix)
      dmat = dispersal_multipliers_matrix
      maxent01s_param = spPmat_inputs$maxent01s_param
      maxent01v_param = spPmat_inputs$maxent01v_param
      maxent01j_param = spPmat_inputs$maxent01j_param
      maxent01y_param = spPmat_inputs$maxent01y_param
      l = spPmat_inputs$l
      numareas = max(sapply(X = spPmat_inputs$l, FUN = length), 
                     na.rm = TRUE) + 0
      maxent01s = relative_probabilities_of_subsets(max_numareas = numareas, 
                                                    maxent_constraint_01 = maxent01s_param, NA_val = 0)
      maxent01v = relative_probabilities_of_vicariants(max_numareas = numareas, 
                                                       maxent_constraint_01v = maxent01v_param, NA_val = 0)
      maxent01j = relative_probabilities_of_subsets(max_numareas = numareas, 
                                                    maxent_constraint_01 = maxent01j_param, NA_val = 0)
      maxent01y = relative_probabilities_of_subsets(max_numareas = numareas, 
                                                    maxent_constraint_01 = maxent01y_param, NA_val = 0)
      maxprob_as_function_of_ancsize_and_decsize = mapply(FUN = max, 
                                                          maxent01s, maxent01v, maxent01j, maxent01y, 
                                                          MoreArgs = list(na.rm = TRUE))
      maxprob_as_function_of_ancsize_and_decsize = matrix(data = maxprob_as_function_of_ancsize_and_decsize, 
                                                          nrow = nrow(maxent01s), ncol = ncol(maxent01s))
      maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize > 
                                                   0] = 1
      maxprob_as_function_of_ancsize_and_decsize[maxprob_as_function_of_ancsize_and_decsize <= 
                                                   0] = 0
      max_minsize_as_function_of_ancsize = apply(X = maxprob_as_function_of_ancsize_and_decsize, 
                                                 MARGIN = 1, FUN = maxsize)
      if (include_null_range == TRUE) {
        state_space_size_Qmat_to_cladoMat = -1
      }
      else {
        state_space_size_Qmat_to_cladoMat = 0
      }
      tmpca_1 = rep(1, sum(states_allowed_TF) + state_space_size_Qmat_to_cladoMat)
      tmpcb_1 = rep(1, sum(states_allowed_TF) + state_space_size_Qmat_to_cladoMat)
      COO_weights_columnar = rcpp_calc_anclikes_sp_COOweights_faster(Rcpp_leftprobs = tmpca_1, 
                                                                     Rcpp_rightprobs = tmpcb_1, l = l, s = spPmat_inputs$s, 
                                                                     v = spPmat_inputs$v, j = spPmat_inputs$j, y = spPmat_inputs$y, 
                                                                     dmat = dmat, maxent01s = maxent01s, maxent01v = maxent01v, 
                                                                     maxent01j = maxent01j, maxent01y = maxent01y, 
                                                                     max_minsize_as_function_of_ancsize = max_minsize_as_function_of_ancsize, 
                                                                     printmat = FALSE, m = m)
      Rsp_rowsums = rcpp_calc_rowsums_for_COOweights_columnar(COO_weights_columnar = COO_weights_columnar)
      cppSpMethod = 3
      if (exists("COO_weights_columnar") == FALSE) {
        stop("\nERROR_A: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
      }
      if (exists("Rsp_rowsums") == FALSE) {
        stop("\nERROR_B: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
      }
      if (cppSpMethod != 3) {
        stop("\nERROR_C: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
      }
      chainsaw_result = inputs$tree_sections_list[[i]]
      inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips = list()
      for (jj in 1:length(chainsaw_result$return_pieces_list)) {
        treepiece = chainsaw_result$return_pieces_list[[jj]]
        if (is.numeric(treepiece)) {
          do_exponentiation = TRUE
          if (i == num_timeperiods) {
            errortxt = "ERROR: In stratified analysis, your tree must start with a root node, not a branch below the root node."
            stop(errortxt)
          }
          subbranch_length = treepiece
          TF1 = inputs$master_table$stratum == i
          TF2 = inputs$master_table$piecenum == jj
          TF3 = inputs$master_table$piececlass == "subbranch"
          TF = (TF1 + TF2 + TF3) == 3
          anc_node_original_tree = inputs$master_table$node[TF]
          rownum = (1:nrow(inputs$master_table))[TF]
          tmp_master_table_row = inputs$master_table[rownum, 
          ]
          if (nrow(tmp_master_table_row) != 1) {
            stoptxt = paste("\n\nFATAL ERROR in stratified loglike UPPASS calculation at i=", 
                            i, "; jj=", jj, "; ", "inputs$master_table$piececlass == \"subbranch\"", 
                            "\nnrow(tmp_master_table_row) should =1 but instead =", 
                            nrow(tmp_master_table_row), "\n", sep = "")
            stop(stoptxt)
          }
          master_tip_time_bp = tmp_master_table_row$time_bp
          time_top = tmp_master_table_row$time_top
          time_bot = tmp_master_table_row$time_bot
          is_fossil = tmp_master_table_row$fossils
          if (master_tip_time_bp < time_top) {
            do_exponentiation = FALSE
          }
          if ((master_tip_time_bp >= time_top) && (master_tip_time_bp < 
                                                   time_bot) && (is_fossil == TRUE)) {
            amount_to_shorten_by = master_tip_time_bp - 
              time_top
            subbranch_length = subbranch_length - amount_to_shorten_by
            do_exponentiation = TRUE
          }
          if (tmp_master_table_row$edge.length < min_branchlength) {
            do_exponentiation = FALSE
          }
          previous_stratum = i + 1
          previous_stratum_TF = inputs$master_table$stratum == 
            previous_stratum
          node_TF = inputs$master_table$node == anc_node_original_tree
          TF = (previous_stratum_TF + node_TF) == 2
          master_table_row_corresponding_to_anctip = inputs$master_table[TF, 
          ]
          previous_treepiece_num = master_table_row_corresponding_to_anctip$piecenum
          previous_treepiece = inputs$tree_sections_list[[previous_stratum]]$return_pieces_list[[previous_treepiece_num]]
          relprobs_at_tips_of_anc_treepiece = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_tips[[previous_treepiece_num]]
          relprobs_at_branch_bottoms_below_tips_from_previous_stratum = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_bottoms_below_tips[[previous_treepiece_num]]
          if (is.numeric(previous_treepiece) == TRUE) {
            ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece
            ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum
          }
          else {
            tipnum_in_previous_treepiece = master_table_row_corresponding_to_anctip$SUBnode
            ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece[tipnum_in_previous_treepiece, 
            ]
            ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum[tipnum_in_previous_treepiece, 
            ]
          }
          if (do_exponentiation == TRUE) {
            if (sparse == FALSE) {
              if (traitTF == FALSE) {
                actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[states_allowed_TF], 
                                                                                              branch_length = subbranch_length, 
                                                                                              Qmat_tmp)
              }
              else {
                actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], 
                                                                                              branch_length = subbranch_length, 
                                                                                              Qmat_tmp)
              }
              if (include_null_range == TRUE) {
                actual_probs_after_forward_exponentiation[1] = 0
              }
              actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation/sum(actual_probs_after_forward_exponentiation)
            }
            else {
              if (traitTF == FALSE) {
                actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[states_allowed_TF], 
                                                                                               branch_length = subbranch_length, 
                                                                                               tmpQmat_in_REXPOKIT_coo_fmt = tmpQmat_in_REXPOKIT_coo_fmt, 
                                                                                               coo_n = numstates_geogtrait, anorm = NULL, 
                                                                                               check_for_0_rows = TRUE)
              }
              else {
                actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], 
                                                                                               branch_length = subbranch_length, 
                                                                                               tmpQmat_in_REXPOKIT_coo_fmt = tmpQmat_in_REXPOKIT_coo_fmt, 
                                                                                               coo_n = numstates_geogtrait, anorm = NULL, 
                                                                                               check_for_0_rows = TRUE)
              }
              if (include_null_range == TRUE) {
                actual_probs_after_forward_exponentiation[1] = 0
              }
              actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation/sum(actual_probs_after_forward_exponentiation)
            }
            if (traitTF == FALSE) {
              actual_probs_after_forward_exponentiation_new = rep(0, 
                                                                  length(states_allowed_TF))
              actual_probs_after_forward_exponentiation_new[states_allowed_TF] = actual_probs_after_forward_exponentiation
            }
            else {
              actual_probs_after_forward_exponentiation_new = rep(0, 
                                                                  length(wTrait_states_allowed_TF))
              actual_probs_after_forward_exponentiation_new[wTrait_states_allowed_TF] = actual_probs_after_forward_exponentiation
            }
            actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation_new
            if (!is.null(states_allowed_TF)) {
              actual_probs_after_forward_exponentiation[states_allowed_TF == 
                                                          FALSE] = 0
              actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation/sum(actual_probs_after_forward_exponentiation)
            }
          }
          else {
            actual_probs_after_forward_exponentiation = ancprobs_at_subbranch_bottom
          }
          if (any(is.na(actual_probs_after_forward_exponentiation))) {
            print("i, jj, anc")
            print(i)
            print(jj)
            print(anc)
            print("actual_probs_after_forward_exponentiation")
            print(actual_probs_after_forward_exponentiation)
            print("ancprobs_at_subbranch_bottom")
            print(ancprobs_at_subbranch_bottom)
            print("ancprobs_at_bottom_of_total_branch")
            print(ancprobs_at_bottom_of_total_branch)
            stop("ERROR #1: see stratified code")
          }
          relprobs_at_tips_for_next_stratum_up = actual_probs_after_forward_exponentiation
          relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = ancprobs_at_bottom_of_total_branch
          inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips[[jj]] = relprobs_at_tips_for_next_stratum_up
          relprobs_at_tips_for_next_stratum_up
          inputs$tree_sections_list[[i]]$pieces_relprobs_at_bottoms_below_tips[[jj]] = relprobs_at_branch_bottoms_below_tips_for_next_stratum_up
        }
        else {
          tmp_subtree = treepiece
          use_fixnodes_on_uppass = TRUE
          if (use_fixnodes_on_uppass) {
            tmp_fixnode = NULL
            tmp_fixlikes = NULL
            if ((!is.null(fixnode)) && (length(fixnode) > 
                                        0)) {
              if (length(fixnode) > 1) {
                TF1 = inputs$master_table$stratum == 
                  i
                TF2 = inputs$master_table$piecenum == 
                  jj
                TF3 = inputs$master_table$piececlass == 
                  "subtree"
                TF = ((TF1 + TF2 + TF3) == 3)
                tmprows = inputs$master_table[TF, ]
                fixnodes_in_subtree_TF = fixnode %in% 
                  tmprows$node
                if (sum(fixnodes_in_subtree_TF) > 0) {
                  temporary_fixnodes = fixnode[fixnodes_in_subtree_TF]
                  temporary_fixlikes = fixlikes[fixnodes_in_subtree_TF, 
                  ]
                  subtree_rows_in_fixnodes_TF = tmprows$node %in% 
                    fixnode
                  subtree_fixnode_master_nodenums = tmprows$node[subtree_rows_in_fixnodes_TF]
                  subtree_fixnode_nums = tmprows$SUBnode[subtree_rows_in_fixnodes_TF]
                  order_subtree_fixnode_nums = order(subtree_fixnode_nums)
                  subtree_fixnode_nums = subtree_fixnode_nums[order_subtree_fixnode_nums]
                  if (length(order_subtree_fixnode_nums) > 
                      1) {
                    temporary_fixlikes = temporary_fixlikes[order_subtree_fixnode_nums, 
                    ]
                  }
                }
                else {
                  temporary_fixnodes = NULL
                  subtree_fixnode_master_nodenums = NULL
                  subtree_fixnode_nums = NULL
                  temporary_fixlikes = NULL
                }
                TF1 = unique(tmprows$stratum) == i
                TF2 = unique(tmprows$piecenum) == jj
                TF3 = unique(tmprows$piececlass) == 
                  "subtree"
                TF = ((TF1 + TF2 + TF3) == 3)
                if (TF == TRUE) {
                  if (length(subtree_fixnode_nums) == 
                      0) {
                    subtree_fixnode_nums = NULL
                    temporary_fixlikes = NULL
                  }
                  tmp_fixnode = subtree_fixnode_nums
                  tmp_fixlikes = temporary_fixlikes
                }
                else {
                  tmp_fixnode = NULL
                  tmp_fixlikes = NULL
                }
              }
              else {
                temporary_fixnode = fixnode
                temporary_fixlikes = c(fixlikes)
                TF1 = inputs$master_table$node == temporary_fixnode
                TF2 = inputs$master_table$SUBnode.type == 
                  "root"
                TF3 = inputs$master_table$SUBnode.type == 
                  "internal"
                TF = ((TF1 + TF2 + TF3) == 2)
                tmprow = inputs$master_table[TF, ]
                TF1 = tmprow$stratum == i
                TF2 = tmprow$piecenum == jj
                TF3 = tmprow$piececlass == "subtree"
                TF = ((TF1 + TF2 + TF3) == 3)
                if (TF == TRUE) {
                  tmp_fixnode = tmprow$SUBnode
                  tmp_fixlikes = temporary_fixlikes
                }
                else {
                  tmp_fixnode = NULL
                  tmp_fixlikes = NULL
                }
              }
            }
          }
          tmp_subtree_tipnums = 1:length(tmp_subtree$tip.label)
          for (iter in 1:length(tmp_subtree_tipnums)) {
            subtree_tip = tmp_subtree_tipnums[iter]
            TF1 = inputs$master_table$stratum == i
            TF2 = inputs$master_table$piecenum == jj
            TF3 = inputs$master_table$piececlass == 
              "subtree"
            TF4 = inputs$master_table$SUBnode == subtree_tip
            TF = (TF1 + TF2 + TF3 + TF4) == 4
            rownum = (1:nrow(inputs$master_table))[TF]
            tmp_master_table_row = inputs$master_table[rownum, 
            ]
            if (nrow(tmp_master_table_row) != 1) {
              stoptxt = paste("\n\nFATAL ERROR in stratified loglike UPPASS calculation at i=", 
                              i, "; jj=", jj, "; ", "inputs$master_table$piececlass == \"subtree\"", 
                              "; subtree_tip=", subtree_tip, "\nnrow(tmp_master_table_row) should =1 but instead =", 
                              nrow(tmp_master_table_row), "\n", sep = "")
              stop(stoptxt)
            }
            master_tip_time_bp = tmp_master_table_row$time_bp
            time_top = tmp_master_table_row$time_top
            time_bot = tmp_master_table_row$time_bot
            is_fossil = tmp_master_table_row$fossils
            if ((master_tip_time_bp >= time_top) && 
                (master_tip_time_bp < time_bot) && is_fossil == 
                TRUE) {
              amount_to_shorten_by = master_tip_time_bp - 
                time_top
              tmp2_edgeTF = tmp_subtree$edge[, 2] == 
                subtree_tip
              tmp2_edgenum = (1:nrow(tmp_subtree$edge))[tmp2_edgeTF]
            }
          }
          phy2 <- reorder(tmp_subtree, "pruningwise")
          tipnames = phy2$tip.label
          num_internal_nodes = phy2$Nnode
          edges_to_visit_uppass = seq(from = (num_internal_nodes * 
                                                2), by = -2, length.out = num_internal_nodes)
          tmpj = edges_to_visit_uppass[1]
          tmpi = tmpj - 1
          anc <- phy2$edge[tmpi, 1]
          numnodes = num_internal_nodes + length(tmp_subtree$tip.label)
          TFi = inputs$master_table$stratum == i
          TFjj = inputs$master_table$piecenum == jj
          TF_SUBnode = inputs$master_table$SUBnode == 
            anc
          TF = ((TFi + TFjj + TF_SUBnode) == 3)
          anc_node_original_tree = inputs$master_table$node[TF]
          if (traitTF == FALSE) {
            tmp_relprobs_at_branchtop_AT_node_UPPASS = matrix(data = 0, 
                                                              nrow = numnodes, sum(states_allowed_TF))
            tmp_relprobs_at_branchbot_BELOW_node_UPPASS = matrix(data = 0, 
                                                                 nrow = numnodes, sum(states_allowed_TF))
          }
          else {
            tmp_relprobs_at_branchtop_AT_node_UPPASS = matrix(data = 0, 
                                                              nrow = numnodes, sum(states_allowed_TF) * 
                                                                num_trait_states)
            tmp_relprobs_at_branchbot_BELOW_node_UPPASS = matrix(data = 0, 
                                                                 nrow = numnodes, sum(states_allowed_TF) * 
                                                                   num_trait_states)
          }
          master_tree_nodenums = NULL
          for (rownum in 1:nrow(tmp_relprobs_at_branchtop_AT_node_UPPASS)) {
            tmp_relprobs = tmp_relprobs_at_branchtop_AT_node_UPPASS[rownum, 
            ]
            subtree_node = rownum
            TF1 = inputs$master_table$stratum == i
            TF2 = inputs$master_table$piecenum == jj
            TF3 = inputs$master_table$piececlass == 
              "subtree"
            TF4 = inputs$master_table$SUBnode == subtree_node
            TF5a = inputs$master_table$node.type == 
              "internal"
            TF5b = inputs$master_table$node.type == 
              "tip"
            TF5c = inputs$master_table$node.type == 
              "root"
            TF_subtrees = (TF1 + TF2 + TF3 + TF4 + TF5a + 
                             TF5b + TF5c) == 5
            master_tree_nodenums = c(master_tree_nodenums, 
                                     inputs$master_table$node[TF_subtrees])
          }
          if (traitTF == FALSE) {
            tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[master_tree_nodenums, 
            ][, states_allowed_TF]
          }
          else {
            tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[master_tree_nodenums, 
            ][, wTrait_states_allowed_TF]
          }
          if (sum(states_allowed_TF) == 1) {
            if (traitTF == FALSE) {
              tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = matrix(data = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS, 
                                                                     ncol = 1)
            }
            else {
              tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS = matrix(data = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS, 
                                                                     ncol = sum(wTrait_states_allowed_TF))
            }
          }
          if (i == num_timeperiods) {
            if (traitTF == FALSE) {
              tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
              ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, 
              ][states_allowed_TF]
            }
            else {
              tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
              ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, 
              ][wTrait_states_allowed_TF]
            }
          }
          else {
            if ((is.numeric(phy2$root.edge) == TRUE) && 
                (!is.null(phy2$root.edge)) && (phy2$root.edge > 
                                               0)) {
              root_edge_length = phy2$root.edge
              previous_stratum = i + 1
              previous_stratum_TF = inputs$master_table$stratum == 
                previous_stratum
              node_TF = inputs$master_table$node == 
                anc_node_original_tree
              TF = (previous_stratum_TF + node_TF) == 
                2
              master_table_row_corresponding_to_anctip = inputs$master_table[TF, 
              ]
              previous_treepiece_num = master_table_row_corresponding_to_anctip$piecenum
              previous_treepiece = inputs$tree_sections_list[[previous_stratum]]$return_pieces_list[[previous_treepiece_num]]
              relprobs_at_tips_of_anc_treepiece = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_tips[[previous_treepiece_num]]
              relprobs_at_branch_bottoms_below_tips_from_previous_stratum = inputs$tree_sections_list[[previous_stratum]]$pieces_relprobs_at_bottoms_below_tips[[previous_treepiece_num]]
              if (is.numeric(previous_treepiece) == 
                  TRUE) {
                ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece
                ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum
              }
              else {
                tipnum_in_previous_treepiece = master_table_row_corresponding_to_anctip$SUBnode
                ancprobs_at_subbranch_bottom = relprobs_at_tips_of_anc_treepiece[tipnum_in_previous_treepiece, 
                ]
                ancprobs_at_bottom_of_total_branch = relprobs_at_branch_bottoms_below_tips_from_previous_stratum[tipnum_in_previous_treepiece, 
                ]
              }
              if (sparse == FALSE) {
                if (traitTF == FALSE) {
                  actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[states_allowed_TF], 
                                                                                                branch_length = root_edge_length, 
                                                                                                Qmat_tmp)
                }
                else {
                  actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_dense(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], 
                                                                                                branch_length = root_edge_length, 
                                                                                                Qmat_tmp)
                }
              }
              else {
                if (traitTF == FALSE) {
                  actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[states_allowed_TF], 
                                                                                                 branch_length = root_edge_length, 
                                                                                                 tmpQmat_in_REXPOKIT_coo_fmt = tmpQmat_in_REXPOKIT_coo_fmt, 
                                                                                                 coo_n = numstates_geogtrait, anorm = NULL, 
                                                                                                 check_for_0_rows = TRUE)
                }
                else {
                  actual_probs_after_forward_exponentiation = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom = ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF], 
                                                                                                 branch_length = root_edge_length, 
                                                                                                 tmpQmat_in_REXPOKIT_coo_fmt = tmpQmat_in_REXPOKIT_coo_fmt, 
                                                                                                 coo_n = numstates_geogtrait, anorm = NULL, 
                                                                                                 check_for_0_rows = TRUE)
                }
              }
              if (traitTF == FALSE) {
                actual_probs_after_forward_exponentiation_new = rep(0, 
                                                                    length(states_allowed_TF))
                actual_probs_after_forward_exponentiation_new[states_allowed_TF] = actual_probs_after_forward_exponentiation
                actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation_new
              }
              else {
                actual_probs_after_forward_exponentiation_new = rep(0, 
                                                                    length(wTrait_states_allowed_TF))
                actual_probs_after_forward_exponentiation_new[wTrait_states_allowed_TF] = actual_probs_after_forward_exponentiation
                actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation_new
              }
              if (include_null_range == TRUE) {
                actual_probs_after_forward_exponentiation[1] = 0
              }
              actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation/sum(actual_probs_after_forward_exponentiation)
              if (!is.null(states_allowed_TF)) {
                actual_probs_after_forward_exponentiation[states_allowed_TF == 
                                                            FALSE] = 0
                actual_probs_after_forward_exponentiation = actual_probs_after_forward_exponentiation/sum(actual_probs_after_forward_exponentiation)
              }
              ancprobs_at_subtree_root = actual_probs_after_forward_exponentiation
              if (traitTF == FALSE) {
                tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                ] = actual_probs_after_forward_exponentiation[states_allowed_TF]
                tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
                ] = ancprobs_at_bottom_of_total_branch[states_allowed_TF]
              }
              else {
                tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                ] = actual_probs_after_forward_exponentiation[wTrait_states_allowed_TF]
                tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
                ] = ancprobs_at_bottom_of_total_branch[wTrait_states_allowed_TF]
              }
            }
            else {
              if (traitTF == FALSE) {
                tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, 
                ][states_allowed_TF]
                tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
                ] = ancprobs_at_bottom_of_total_branch[states_allowed_TF]
              }
              else {
                tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                ] = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[anc_node_original_tree, 
                ][wTrait_states_allowed_TF]
                tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
                ] = ancprobs_at_bottom_of_total_branch[wTrait_states_allowed_TF]
              }
            }
            if (any(is.na(tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
            ]))) {
              print("i, jj, anc")
              print(i)
              print(jj)
              print(anc)
              print("actual_probs_after_forward_exponentiation")
              print(actual_probs_after_forward_exponentiation)
              print("tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, ]")
              print(tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
              ])
              print("anc_node_original_tree")
              print(anc_node_original_tree)
              print("states_allowed_TF")
              print(states_allowed_TF)
              if (traitTF == FALSE) {
                print("ancprobs_at_subbranch_bottom[states_allowed_TF]")
                print(ancprobs_at_subbranch_bottom[states_allowed_TF])
              }
              else {
                print("wTrait_states_allowed_TF")
                print(wTrait_states_allowed_TF)
                print("ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF]")
                print(ancprobs_at_subbranch_bottom[wTrait_states_allowed_TF])
              }
              stop("ERROR #2: see stratified code")
            }
            if (any(is.na(tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
            ]))) {
              print("i, jj, anc")
              print(i)
              print(jj)
              print(anc)
              print("actual_probs_after_forward_exponentiation")
              print(actual_probs_after_forward_exponentiation)
              print("ancprobs_at_bottom_of_total_branch")
              print(ancprobs_at_bottom_of_total_branch)
              print("tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, ]")
              print(tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
              ])
              stop("ERROR #3: see stratified code")
            }
          }
          if (sparse == FALSE) {
            independent_likelihoods_on_each_branch = vector("list", 
                                                            length(phy2$edge.length))
            tmpmatrix = matrix(data = 0, nrow = nrow(Qmat_tmp), 
                               ncol = ncol(Qmat_tmp))
            for (m in 1:length(phy2$edge.length)) {
              independent_likelihoods_on_each_branch[[m]] = tmpmatrix
            }
            if (!is.null(cluster_already_open)) {
              if (.Platform$GUI == "AQUA") {
                cat("In calc_loglike_sp(), cluster_already_open=", 
                    cluster_already_open, " which means you want to calculate likelihoods on branches using a multicore option.\n", 
                    sep = "")
                cat("But .Platform$GUI='AQUA', which means you are running the Mac GUI R.app version of R.  Parallel multicore functions, e.g. as accessed via \n", 
                    sep = "")
                cat("library(parallel), are apparently dangerous/will crash R.app (google multicore 'R.app').  So, changing to cluster_already_open=NULL.\n", 
                    sep = "")
                cluster_already_open = NULL
              }
            }
            if (!is.null(cluster_already_open)) {
              independent_likelihoods_on_each_branch = clusterApply(cl = cluster_already_open, 
                                                                    x = phy2$edge.length, fun = expokit_dgpadm_Qmat2, 
                                                                    Qmat = Qmat_tmp, transpose_needed = TRUE)
            }
            else {
              independent_likelihoods_on_each_branch = mapply_likelihoods(Qmat_tmp, 
                                                                          phy2, transpose_needed = TRUE)
            }
          }
          rootnode = length(phy2$tip.label) + 1
          for (uj in edges_to_visit_uppass) {
            ui <- uj - 1
            left_desc_nodenum <- phy2$edge[ui, 2]
            right_desc_nodenum <- phy2$edge[uj, 2]
            anc <- phy2$edge[ui, 1]
            anc_edgenum_TF = phy2$edge[, 2] == anc
            anc_edgenum = (1:length(anc_edgenum_TF))[anc_edgenum_TF]
            mother_of_anc_TF = phy2$edge[, 2] == anc
            mother_of_anc = phy2$edge[mother_of_anc_TF, 
                                      1]
            sister_of_anc_TF = phy2$edge[, 1] == mother_of_anc
            sister_of_anc_TF2 = (sister_of_anc_TF + 
                                   mother_of_anc_TF) == 1
            sister_of_anc = phy2$edge[sister_of_anc_TF2, 
                                      2]
            mother_of_anc
            sister_of_anc
            sister_is_LR = "rootnode"
            if (anc != rootnode) {
              if (sister_of_anc > anc) {
                sister_is_LR = "right"
              }
              else {
                sister_is_LR = "left"
              }
            }
            left_edge_TF = phy2$edge[, 2] == left_desc_nodenum
            right_edge_TF = phy2$edge[, 2] == right_desc_nodenum
            left_edgenum = (1:length(left_edge_TF))[left_edge_TF]
            right_edgenum = (1:length(right_edge_TF))[right_edge_TF]
            is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < 
              min_branchlength
            is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < 
              min_branchlength
            hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 
              0
            tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
            ]
            sum_states_allowed = sum(states_allowed_TF)
            if (traitTF == TRUE) {
              wTrait_sum_states_allowed = sum(states_allowed_TF) * 
                num_trait_states
            }
            if (hooknode_TF == TRUE) {
              temp_COO_weights_columnar = COO_weights_columnar
              if (include_null_range == TRUE) {
                highest_clado_state_0based_considering_null_range = sum_states_allowed - 
                  2
              }
              else {
                highest_clado_state_0based_considering_null_range = sum_states_allowed - 
                  1
              }
              temp_COO_weights_columnar[[1]] = 0:highest_clado_state_0based_considering_null_range
              temp_COO_weights_columnar[[2]] = 0:highest_clado_state_0based_considering_null_range
              temp_COO_weights_columnar[[3]] = 0:highest_clado_state_0based_considering_null_range
              temp_COO_weights_columnar[[4]] = rep(1, 
                                                   highest_clado_state_0based_considering_null_range + 
                                                     1)
            }
            else {
              temp_COO_weights_columnar = COO_weights_columnar
            }
            num_nonzero_split_scenarios = length(COO_weights_columnar[[1]])
            TFi = inputs$master_table$stratum == i
            TFjj = inputs$master_table$piecenum == jj
            TF_SUBnode = inputs$master_table$SUBnode == 
              anc
            TF = ((TFi + TFjj + TF_SUBnode) == 3)
            anc_node_original_tree = inputs$master_table$node[TF]
            global_root_TF = inputs$master_table$node.type[TF]
            if ((anc == rootnode) && (global_root_TF == 
                                      TRUE)) {
              probs_at_mother = 1/length(starting_probs)
              likes_at_sister = 1/length(starting_probs)
              left_branch_downpass_likes = NULL
              right_branch_downpass_likes = NULL
              probs_of_mother_and_sister_uppass_to_anc = 1/length(starting_probs)
              tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
              ] = probs_of_mother_and_sister_uppass_to_anc
            }
            else {
              if (anc == rootnode) {
                probs_at_mother = 1/length(starting_probs)
                probs_of_mother_and_sister_uppass_to_anc = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                ]
              }
              else {
                probs_at_mother = tmp_relprobs_at_branchtop_AT_node_UPPASS[mother_of_anc, 
                ]
              }
              likes_at_sister_branch_bottom = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS[sister_of_anc, 
              ]
              if (sister_is_LR == "left") {
                left_branch_downpass_likes = likes_at_sister_branch_bottom
                right_branch_downpass_likes = NULL
              }
              if (sister_is_LR == "right") {
                left_branch_downpass_likes = NULL
                right_branch_downpass_likes = likes_at_sister_branch_bottom
              }
              if (anc != rootnode) {
                if (traitTF == FALSE) {
                  uppass_probs_at_bottom_below_anc_results = calc_uppass_probs_new2(probs_ancstate = probs_at_mother, 
                                                                                    COO_weights_columnar = temp_COO_weights_columnar, 
                                                                                    numstates = sum_states_allowed, 
                                                                                    include_null_range = include_null_range, 
                                                                                    left_branch_downpass_likes = left_branch_downpass_likes, 
                                                                                    right_branch_downpass_likes = right_branch_downpass_likes, 
                                                                                    Rsp_rowsums = NULL)
                }
                else {
                  uppass_probs_at_bottom_below_anc_results = calc_uppass_probs_new2(probs_ancstate = probs_at_mother, 
                                                                                    COO_weights_columnar = temp_COO_weights_columnar, 
                                                                                    numstates = wTrait_sum_states_allowed, 
                                                                                    include_null_range = include_null_range, 
                                                                                    left_branch_downpass_likes = left_branch_downpass_likes, 
                                                                                    right_branch_downpass_likes = right_branch_downpass_likes, 
                                                                                    Rsp_rowsums = NULL)
                }
                if (sister_is_LR == "left") {
                  Rprobs_brbot_below_anc = uppass_probs_at_bottom_below_anc_results$relprobs_just_after_speciation_UPPASS_Right
                  tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
                  ] = Rprobs_brbot_below_anc
                }
                if (sister_is_LR == "right") {
                  Lprobs_brbot_below_anc = uppass_probs_at_bottom_below_anc_results$relprobs_just_after_speciation_UPPASS_Left
                  tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
                  ] = Lprobs_brbot_below_anc
                }
                probs_at_branch_bot = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[anc, 
                ]
                if (force_sparse == FALSE) {
                  probs_of_mother_and_sister_uppass_to_anc = probs_at_branch_bot %*% 
                    expokit_dgpadm_Qmat2(times = phy2$edge.length[anc_edgenum], 
                                         Qmat = Qmat_tmp, transpose_needed = TRUE)
                }
                else {
                  probs_of_mother_and_sister_uppass_to_anc = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom = probs_at_branch_bot, 
                                                                                                branch_length = phy2$edge.length[anc_edgenum], 
                                                                                                tmpQmat_in_REXPOKIT_coo_fmt = tmpQmat_in_REXPOKIT_coo_fmt, 
                                                                                                coo_n = length(probs_at_branch_bot), 
                                                                                                anorm = NULL, check_for_0_rows = TRUE)
                }
              }
              else {
                probs_of_mother_and_sister_uppass_to_anc
              }
              tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
              ] = probs_of_mother_and_sister_uppass_to_anc
            }
            if (left_desc_nodenum <= length(phy2$tip.label)) {
              probs_at_anc = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
              ]
              left_branch_downpass_likes = NULL
              right_branch_downpass_likes = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS[right_desc_nodenum, 
              ]
              if (traitTF == FALSE) {
                uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate = probs_at_anc, 
                                                                                  COO_weights_columnar = temp_COO_weights_columnar, 
                                                                                  numstates = sum_states_allowed, include_null_range = include_null_range, 
                                                                                  left_branch_downpass_likes = left_branch_downpass_likes, 
                                                                                  right_branch_downpass_likes = right_branch_downpass_likes, 
                                                                                  Rsp_rowsums = NULL)
              }
              else {
                uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate = probs_at_anc, 
                                                                                  COO_weights_columnar = temp_COO_weights_columnar, 
                                                                                  numstates = wTrait_sum_states_allowed, 
                                                                                  include_null_range = include_null_range, 
                                                                                  left_branch_downpass_likes = left_branch_downpass_likes, 
                                                                                  right_branch_downpass_likes = right_branch_downpass_likes, 
                                                                                  Rsp_rowsums = NULL)
              }
              Lprobs_brbot_below_tip = uppass_probs_at_bottom_below_tip_results$relprobs_just_after_speciation_UPPASS_Left
              if (force_sparse == FALSE) {
                Lprobs_brtop_AT_tip = Lprobs_brbot_below_tip %*% 
                  expokit_dgpadm_Qmat2(times = phy2$edge.length[left_edgenum], 
                                       Qmat = Qmat_tmp, transpose_needed = TRUE)
              }
              else {
                Lprobs_brtop_AT_tip = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom = Lprobs_brbot_below_tip, 
                                                                         branch_length = phy2$edge.length[left_edgenum], 
                                                                         tmpQmat_in_REXPOKIT_coo_fmt = tmpQmat_in_REXPOKIT_coo_fmt, 
                                                                         coo_n = numstates_geogtrait, anorm = NULL, 
                                                                         check_for_0_rows = TRUE)
              }
              tmp_relprobs_at_branchbot_BELOW_node_UPPASS[left_desc_nodenum, 
              ] = Lprobs_brbot_below_tip
              tmp_relprobs_at_branchtop_AT_node_UPPASS[left_desc_nodenum, 
              ] = Lprobs_brtop_AT_tip
            }
            if (right_desc_nodenum <= length(phy2$tip.label)) {
              probs_at_anc = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
              ]
              right_branch_downpass_likes = NULL
              left_branch_downpass_likes = tmp_relprobs_at_branchbot_BELOW_node_DOWNPASS[left_desc_nodenum, 
              ]
              if (traitTF == FALSE) {
                uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate = probs_at_anc, 
                                                                                  COO_weights_columnar = temp_COO_weights_columnar, 
                                                                                  numstates = sum_states_allowed, include_null_range = include_null_range, 
                                                                                  right_branch_downpass_likes = right_branch_downpass_likes, 
                                                                                  left_branch_downpass_likes = left_branch_downpass_likes, 
                                                                                  Rsp_rowsums = NULL)
              }
              else {
                uppass_probs_at_bottom_below_tip_results = calc_uppass_probs_new2(probs_ancstate = probs_at_anc, 
                                                                                  COO_weights_columnar = temp_COO_weights_columnar, 
                                                                                  numstates = wTrait_sum_states_allowed, 
                                                                                  include_null_range = include_null_range, 
                                                                                  right_branch_downpass_likes = right_branch_downpass_likes, 
                                                                                  left_branch_downpass_likes = left_branch_downpass_likes, 
                                                                                  Rsp_rowsums = NULL)
              }
              Rprobs_brbot_below_tip = uppass_probs_at_bottom_below_tip_results$relprobs_just_after_speciation_UPPASS_Right
              if (force_sparse == FALSE) {
                Rprobs_brtop_AT_tip = Rprobs_brbot_below_tip %*% 
                  expokit_dgpadm_Qmat2(times = phy2$edge.length[right_edgenum], 
                                       Qmat = Qmat_tmp, transpose_needed = TRUE)
              }
              else {
                Rprobs_brtop_AT_tip = calc_prob_forward_onebranch_sparse(relprobs_branch_bottom = Rprobs_brbot_below_tip, 
                                                                         branch_length = phy2$edge.length[right_edgenum], 
                                                                         tmpQmat_in_REXPOKIT_coo_fmt = tmpQmat_in_REXPOKIT_coo_fmt, 
                                                                         coo_n = numstates_geogtrait, anorm = NULL, 
                                                                         check_for_0_rows = TRUE)
              }
              tmp_relprobs_at_branchbot_BELOW_node_UPPASS[right_desc_nodenum, 
              ] = Rprobs_brbot_below_tip
              tmp_relprobs_at_branchtop_AT_node_UPPASS[right_desc_nodenum, 
              ] = Rprobs_brtop_AT_tip
            }
            use_fixnodes_on_uppass = TRUE
            if (use_fixnodes_on_uppass) {
              if (!is.null(fixnode)) {
                if (length(tmp_fixnode) > 1) {
                  TF = (anc == tmp_fixnode)
                  temporary_fixnode = tmp_fixnode[TF]
                  temporary_fixlikes = c(tmp_fixlikes[TF, 
                  ])
                }
                else {
                  temporary_fixnode = tmp_fixnode
                  temporary_fixlikes = c(tmp_fixlikes)
                }
                if ((length(temporary_fixnode) > 0) && 
                    (anc == temporary_fixnode)) {
                  if (traitTF == FALSE) {
                    tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                    ] = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                    ] * temporary_fixlikes[states_allowed_TF]
                  }
                  else {
                    tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                    ] = tmp_relprobs_at_branchtop_AT_node_UPPASS[anc, 
                    ] * temporary_fixlikes[wTrait_states_allowed_TF]
                  }
                }
              }
            }
          }
          for (rownum in 1:nrow(tmp_relprobs_at_branchtop_AT_node_UPPASS)) {
            tmp_relprobs = tmp_relprobs_at_branchtop_AT_node_UPPASS[rownum, 
            ]
            subtree_node = rownum
            TF1 = inputs$master_table$stratum == i
            TF2 = inputs$master_table$piecenum == jj
            TF3 = inputs$master_table$piececlass == 
              "subtree"
            TF4 = inputs$master_table$SUBnode == subtree_node
            TF5a = inputs$master_table$node.type == 
              "internal"
            TF5b = inputs$master_table$node.type == 
              "tip"
            TF5c = inputs$master_table$node.type == 
              "root"
            TF_subtrees = (TF1 + TF2 + TF3 + TF4 + TF5a + 
                             TF5b + TF5c) == 5
            TF = TF_subtrees
            relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum = inputs$master_table$node[TF]
            if (traitTF == FALSE) {
              relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, 
              ][states_allowed_TF] = tmp_relprobs
            }
            else {
              relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, 
              ][wTrait_states_allowed_TF] = tmp_relprobs
            }
            tmp_relprobs = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[rownum, 
            ]
            if (traitTF == FALSE) {
              relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, 
              ][states_allowed_TF] = tmp_relprobs
            }
            else {
              relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS[relative_probs_of_each_state_at_branch_top_AT_node_UPPASS_rownum, 
              ][wTrait_states_allowed_TF] = tmp_relprobs
            }
          }
          if (traitTF == FALSE) {
            relprobs_at_tips_for_next_stratum_up = matrix(0, 
                                                          nrow = length(phy2$tip.label), ncol = length(states_allowed_TF))
            relprobs_at_tips_for_next_stratum_up[, states_allowed_TF] = tmp_relprobs_at_branchtop_AT_node_UPPASS[1:length(phy2$tip.label), 
            ]
            relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = matrix(0, 
                                                                               nrow = length(phy2$tip.label), ncol = length(states_allowed_TF))
            relprobs_at_branch_bottoms_below_tips_for_next_stratum_up[, 
                                                                      states_allowed_TF] = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[1:length(phy2$tip.label), 
                                                                      ]
          }
          else {
            relprobs_at_tips_for_next_stratum_up = matrix(0, 
                                                          nrow = length(phy2$tip.label), ncol = length(wTrait_states_allowed_TF))
            relprobs_at_tips_for_next_stratum_up[, wTrait_states_allowed_TF] = tmp_relprobs_at_branchtop_AT_node_UPPASS[1:length(phy2$tip.label), 
            ]
            relprobs_at_branch_bottoms_below_tips_for_next_stratum_up = matrix(0, 
                                                                               nrow = length(phy2$tip.label), ncol = length(wTrait_states_allowed_TF))
            relprobs_at_branch_bottoms_below_tips_for_next_stratum_up[, 
                                                                      wTrait_states_allowed_TF] = tmp_relprobs_at_branchbot_BELOW_node_UPPASS[1:length(phy2$tip.label), 
                                                                      ]
          }
          inputs$tree_sections_list[[i]]$pieces_relprobs_at_tips[[jj]] = relprobs_at_tips_for_next_stratum_up
          inputs$tree_sections_list[[i]]$pieces_relprobs_at_bottoms_below_tips[[jj]] = relprobs_at_branch_bottoms_below_tips_for_next_stratum_up
        }
      }
    }
    tipnums_of_master_tree = 1:length(original_phy$tip.labe)
    for (tn in 1:length(tipnums_of_master_tree)) {
      TF1 = inputs$master_table$piececlass == "orig_tip"
      TF2 = inputs$master_table$node == tn
      TF = ((TF1 + TF2) == 2)
      tmprow = inputs$master_table[TF, ]
      TF1 = inputs$master_table$node == tmprow$node
      TF2 = inputs$master_table$time_top == tmprow$time_top
      TF3 = inputs$master_table$piececlass != "orig_tip"
      TF = ((TF1 + TF2 + TF3) == 3)
      tmprow2 = inputs$master_table[TF, ]
      tmp_stratum = tmprow2$stratum
      tmp_piecenum = tmprow2$piecenum
      if (tmprow2$piececlass == "subtree") {
        tmp_tipnum = tmprow2$SUBnode
        tmp_tipprobs_at_top_UPPASS = inputs$tree_sections_list[[tmp_stratum]]$pieces_relprobs_at_tips[[tmp_piecenum]][tmp_tipnum, 
        ]
      }
      if (tmprow2$piececlass == "subbranch") {
        tmp_tipnum = tmprow2$SUBnode
        tmp_tipprobs_at_top_UPPASS = inputs$tree_sections_list[[tmp_stratum]]$pieces_relprobs_at_tips[[tmp_piecenum]]
      }
      relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn, 
      ] = tmp_tipprobs_at_top_UPPASS
      if (any(is.na(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn, 
      ]))) {
        print("i, jj, anc, tn")
        print(i)
        print(jj)
        print(anc)
        print(tn)
        print("relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn,]")
        print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS[tn, 
        ])
        print("tmp_tipprobs_at_top_UPPASS")
        print(tmp_tipprobs_at_top_UPPASS)
        print("tmprow2$piececlass")
        print(tmprow2$piececlass)
        stop("ERROR #4: see stratified code")
      }
    }
    cat("\nUppass completed for (STRATIFIED) marginal ancestral states estimation!\n", 
        sep = "")
  }
  if ((return_condlikes_table == TRUE) || (calc_TTL_loglike_from_condlikes_table == 
                                           TRUE)) {
    calc_loglike_sp_stratified_results = NULL
    calc_loglike_sp_stratified_results$final_all_condlikes_of_each_state = final_all_condlikes_of_each_state
    calc_loglike_sp_stratified_results$condlikes_table = condlikes_table
    if (calc_ancprobs == TRUE) {
      calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE
      calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_top_AT_node_UPPASS = relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
      calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS = relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
      ML_marginal_prob_each_state_at_branch_bottom_below_node = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS * 
        relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS
      anc_row_of_master_table_TF = inputs$master_table$node.type == 
        "root"
      anc_node_original_tree = inputs$master_table$node[anc_row_of_master_table_TF]
      anc_node_original_tree
      ML_marginal_prob_each_state_at_branch_bottom_below_node[anc_node_original_tree, 
      ] = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[anc_node_original_tree, 
      ]
      ML_marginal_prob_each_state_at_branch_bottom_below_node
      rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
      ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node/rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
      ML_marginal_prob_each_state_at_branch_bottom_below_node
      sum_MLs_bot = rowSums(ML_marginal_prob_each_state_at_branch_bottom_below_node)
      sum_MLs_bot
      NaN_TF = is.nan(sum_MLs_bot)
      numNaNs = sum(NaN_TF)
      numNaNs
      if (numNaNs > 0) {
        nannodenums = (1:length(NaN_TF))[NaN_TF]
        nannodenums_txt = paste(nannodenums, collapse = ", ", 
                                sep = "")
        txt = paste("\n\nWARNING! ML marginal states at branch bottoms produced ", 
                    numNaNs, " NaNs for nodes:\n", nannodenums_txt, 
                    "\n", "This probably means your downpass probabilities resulted in all 0 probabilities for the node.\n", 
                    "This might occur in a highly constrained model, or if your data strongly contradicts your manual fixed\n", 
                    "likelihoods ('fixlikes') at some node(s) ('fixnode').\n", 
                    "As a 'fix', the downpass probabilities are being used for those nodes. But this is NOT RECOMMENDED!\n", 
                    "You should instead figure out what is causing the problem.", 
                    sep = "")
        cat(txt)
        cat("\n\nPrinting (partial) downpass, uppass, and probability matrices to screen:\n\n", 
            sep = "")
        print(relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS)
        print(relative_probs_of_each_state_at_branch_bottom_below_node_UPPASS)
        print(ML_marginal_prob_each_state_at_branch_bottom_below_node)
        ML_marginal_prob_each_state_at_branch_bottom_below_node[NaN_TF, 
        ] = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[NaN_TF, 
        ]
      }
      ML_marginal_prob_each_state_at_branch_top_AT_node = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS * 
        relative_probs_of_each_state_at_branch_top_AT_node_UPPASS
      ML_marginal_prob_each_state_at_branch_top_AT_node[anc_node_original_tree, 
      ] = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[anc_node_original_tree, 
      ]
      ML_marginal_prob_each_state_at_branch_top_AT_node
      rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
      ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node/rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
      ML_marginal_prob_each_state_at_branch_top_AT_node
      sum_MLs_top = rowSums(ML_marginal_prob_each_state_at_branch_top_AT_node)
      sum_MLs_top
      NaN_TF = is.nan(sum_MLs_top)
      numNaNs = sum(NaN_TF)
      numNaNs
      if (numNaNs > 0) {
        nannodenums = (1:length(NaN_TF))[NaN_TF]
        nannodenums_txt = paste(nannodenums, collapse = ", ", 
                                sep = "")
        txt = paste("\n\nWARNING! ML marginal states at branch tops produced ", 
                    numNaNs, " NaNs for nodes:\n", nannodenums_txt, 
                    "\n", "This probably means your downpass probabilities resulted in all 0 probabilities for the node.\n", 
                    "This might occur in a highly constrained model, or if your data strongly contradicts your manual fixed\n", 
                    "likelihoods ('fixlikes') at some node(s) ('fixnode').\n", 
                    "As a 'fix', the downpass probabilities are being used for those nodes. But this is NOT RECOMMENDED!\n", 
                    "You should instead figure out what is causing the problem.", 
                    sep = "")
        cat(txt)
        cat("\n\nPrinting (partial) downpass, uppass, and probability matrices to screen:\n\n", 
            sep = "")
        print(relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS)
        print(relative_probs_of_each_state_at_branch_top_AT_node_UPPASS)
        print(ML_marginal_prob_each_state_at_branch_top_AT_node)
        ML_marginal_prob_each_state_at_branch_top_AT_node[NaN_TF, 
        ] = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[NaN_TF, 
        ]
      }
      calc_loglike_sp_stratified_results$ML_marginal_prob_each_state_at_branch_bottom_below_node = ML_marginal_prob_each_state_at_branch_bottom_below_node
      calc_loglike_sp_stratified_results$ML_marginal_prob_each_state_at_branch_top_AT_node = ML_marginal_prob_each_state_at_branch_top_AT_node
    }
    calc_loglike_sp_stratified_results$grand_total_likelihood = grand_total_likelihood
    calc_loglike_sp_stratified_results$total_loglikelihood = grand_total_likelihood
    if (calc_ancprobs == TRUE) {
      calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS = relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS
      calc_loglike_sp_stratified_results$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS = relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS
    }
  }
  if (return_what == "loglike") {
    return(grand_total_likelihood)
  }
  if (return_what == "all") {
    return(calc_loglike_sp_stratified_results)
  }
  return(grand_total_likelihood)
  print(paste0("Currently running model ", resfnPLOT, " in ", trfnPH[[1]], " phylogeny.", sep = ""))
}
