# Experimental generic SS3 multinomial length/age composition OSA reweighting.
# Requires: r4ss, afscOSA, SS3 executable, complete runnable model directory.
# Does NOT tune conditional age-at-length, generalized size compositions, or DM parameters.
# OSA-SDNR weighting is a research heuristic, NOT a validated estimator.

.osa_strip <- function(x) trimws(sub('#.*$', '', x))
.osa_starter <- function(path) {
  if (!file.exists(path)) stop('Missing starter.ss: ', path)
  x <- readLines(path, warn=FALSE)
  idx <- which(nzchar(.osa_strip(x)))
  if(length(idx)<2L) stop('Cannot identify data/control filenames in starter.ss')
  vals <- .osa_strip(x[idx[1:2]])
  list(lines=x, index=idx[1:2], data=vals[1], control=vals[2])
}
.osa_rows <- function(path) {
  x <- readLines(path, warn=FALSE)
  head <- grep('^\\s*#\\s*Input variance adjustments factors\\s*:', x, ignore.case=TRUE)
  if(length(head)!=1L) stop('Cannot uniquely locate input variance adjustments block in: ',path)
  tail <- which(grepl('^\\s*-9999\\s+1\\s+0(\\s|$)',x) & seq_along(x)>head)
  if (!length(tail)) stop('Cannot find variance-adjustment terminator in: ',path)
  ind <- seq.int(head+1L,tail[1]-1L)
  ans <- lapply(ind,function(i) {
    s <- .osa_strip(x[i]); if (!nzchar(s)) return(NULL)
    a <- strsplit(s,'[[:space:]]+')[[1]]
    if(length(a)!=3L) return(NULL)
    z <- suppressWarnings(as.numeric(a))
    if(any(!is.finite(z)) || z[1]!=floor(z[1]) || z[2]!=floor(z[2])) return(NULL)
    data.frame(line=i,factor=as.integer(z[1]),fleet=as.integer(z[2]),weight=z[3])
  })
  y <- do.call(rbind,ans)
  if(is.null(y) || !nrow(y)) stop('No variance-adjustment entries found')
  if(anyDuplicated(y[c('factor','fleet')])) stop('Duplicate variance-adjustment entries')
  y
}
.osa_write <- function(src,dest,w) {
  x <- readLines(src,warn=FALSE); rows <- .osa_rows(src)
  for(j in seq_len(nrow(w))) {
    k <- which(rows$factor==w$factor[j] & rows$fleet==w$fleet[j])
    if(length(k)!=1L) stop('Cannot locate control factor/fleet ',w$factor[j], '/',w$fleet[j])
    i <- rows$line[k]; comment <- if(grepl('#',x[i],fixed=TRUE)) paste0(' #',paste(strsplit(x[i],'#',fixed=TRUE)[[1]][-1],collapse='#')) else ''
    x[i] <- sprintf('%d %d %.10g%s',w$factor[j],w$fleet[j],w$weight[j],comment)
  }
  writeLines(x,dest)
}
.osa_residuals <- function(out) {
  z <- out$res
  if(is.matrix(z) && is.numeric(z)) return(as.numeric(z[is.finite(z)]))
  if(!is.data.frame(z)) stop('afscOSA run_osa() returned unexpected $res; inspect its structure')
  c <- intersect(c('resid','osa','OSA','osa_resid','OSA_resid'),names(z))
  if(length(c)!=1L || !is.numeric(z[[c]])) stop('Unexpected afscOSA residual columns: ',paste(names(z),collapse=', '))
  z[[c]][is.finite(z[[c]])]
}
.osa_prepare <- function(tbl, fleet, type, min_N=1, min_years=3) {
  required <- c('Fleet','Yr','Bin','Obs','Exp','Nsamp_adj')
  missing <- setdiff(required,names(tbl))
  if(length(missing)) stop(type, ' composition report missing: ',paste(missing,collapse=', '))
  z <- tbl[!is.na(tbl$Fleet) & tbl$Fleet==fleet,,drop=FALSE]
  if(!nrow(z)) return(list())
  # Only data used in the composition likelihood; fail closed on ambiguous codes.
  if('Use' %in% names(z)) {
    u <- tolower(trimws(as.character(z$Use)))
    z <- z[!is.na(u) & !u %in% c('skip','no','n','0','-1','false'),,drop=FALSE]
  }
  if('Used' %in% names(z)) {
    u <- tolower(trimws(as.character(z$Used)))
    z <- z[!is.na(u) & !u %in% c('skip','no','n','0','-1','false'),,drop=FALSE]
  }
  z <- z[is.finite(suppressWarnings(as.numeric(z$Yr))) &
           is.finite(suppressWarnings(as.numeric(z$Bin))) &
           is.finite(suppressWarnings(as.numeric(z$Nsamp_adj))) &
           as.numeric(z$Nsamp_adj)>=min_N,,drop=FALSE]
  if(!nrow(z)) return(list())
  # Strata that can coexist within a year must NOT be collapsed.
  keys <- intersect(c('Sex','Sexes','Seas','Subseas','Part','Ageerr','Lbin_lo','Lbin_hi',
                       'Area','Method','SupPer','SuprPer'),names(z))
  # A stratum with differing bin structure over years is refused rather than combined.
  tag <- if(length(keys)) do.call(paste,c(lapply(z[keys],as.character),sep='|')) else rep('all',nrow(z))
  groups <- split(z,tag)
  ans <- list()
  for(g in names(groups)) {
    a <- groups[[g]]; a <- a[order(as.numeric(a$Yr),as.numeric(a$Bin)),,drop=FALSE]
    yrs <- sort(unique(as.numeric(a$Yr))); bins <- sort(unique(as.numeric(a$Bin)))
    if(anyDuplicated(a[c('Yr','Bin')]) || nrow(a)!=length(yrs)*length(bins))
      stop('Nonrectangular or duplicated ',type,' compositions, fleet ',fleet,', stratum ',g,
           '. Investigate missing bins/years or additional strata; no automatic pooling.')
    obs <- matrix(as.numeric(a$Obs),nrow=length(yrs),byrow=TRUE)
    exp <- matrix(as.numeric(a$Exp),nrow=length(yrs),byrow=TRUE)
    N <- vapply(yrs,function(y) {
      v <- unique(as.numeric(a$Nsamp_adj[as.numeric(a$Yr)==y]))
      if(length(v)!=1L) stop('Nsamp_adj varies within year in ',type,' fleet ',fleet)
      v
    },numeric(1))
    if(length(yrs)<min_years) {
      warning('Skipping short stratum ',type,' fleet ',fleet, ' [',g,'] (',length(yrs),' years)')
      next
    }
    if(any(!is.finite(obs)) || any(!is.finite(exp)) || any(obs<0) || any(exp<0))
      stop('Invalid composition proportions in ',type,' fleet ',fleet, ' [',g,']')
    so <- rowSums(obs); se <- rowSums(exp)
    if(any(abs(so-1)>0.03) || any(abs(se-1)>0.03))
      stop('Obs/Exp compositions do not sum to 1 in ',type,' fleet ',fleet, ' [',g,']')
    if(any(se<=0) || any(so<=0)) stop('Zero-total composition row')
    ans[[length(ans)+1L]] <- list(obs=obs/so,exp=exp/se,N=N,
                                 index=bins,years=yrs,stratum=g)
  }
  ans
}
.osa_copy_inputs <- function(src,dest,output_root) {
  a <- list.files(src,all.files=TRUE,full.names=TRUE,no..=TRUE)
  a <- a[normalizePath(a,mustWork=FALSE)!=normalizePath(output_root,mustWork=FALSE)]
  # Model root is copied; exclude large / stale SS outputs and prior OSA runs.
  skip <- grepl('^(Report\\.sso|CompReport\\.sso|warning\\.sso|\\.git|OSA_reweight_runs)',basename(a),ignore.case=TRUE)
  a <- a[!skip]
  if(length(a) && !all(file.copy(a,dest,recursive=TRUE,overwrite=TRUE)))
    stop('Could not copy complete model inputs to ',dest)
}

# Combine afscOSA's four native ggplot panels into one vertically stacked PDF page.
# The underlying residuals, data, and individual diagnostic calculations are unchanged.
.osa_composite_page <- function(diagnostic, label, destination, show=FALSE) {
  panel_names <- c('aggcomp', 'qq', 'bubble', 'bubble_pearson')
  missing <- setdiff(panel_names, names(diagnostic))
  if (length(missing)) stop('afscOSA plot_osa() did not provide panels: ',
                            paste(missing, collapse=', '))
  panels <- diagnostic[panel_names]
  if (!all(vapply(panels, inherits, logical(1), 'ggplot')))
    stop('Expected four ggplot panels from afscOSA::plot_osa(). Check your afscOSA version.')
  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(
      nrow=5L, ncol=1L,
      heights=grid::unit(c(0.28, 1.05, 1.20, 1, 1), 'null'))))
    grid::grid.text(label, vp=grid::viewport(layout.pos.row=1L),
                    gp=grid::gpar(fontface='bold', cex=1.0))
    for (i in seq_along(panels)) {
      print(panels[[i]], newpage=FALSE,
            vp=grid::viewport(layout.pos.row=i+1L))
    }
    grid::popViewport()
  }
  grDevices::pdf(destination, width=9, height=15, onefile=FALSE)
  tryCatch(draw(), finally=grDevices::dev.off())
  if (isTRUE(show) && interactive()) draw()
  invisible(destination)
}

# Render the unmodified afscOSA diagnostic components from saved run_osa() objects.
# Each file corresponds to the FINAL FITTED iteration, never proposed next weights.
.osa_plot_final <- function(run_dir, targets, show=FALSE) {
  plot_dir <- file.path(run_dir, 'final_osa_plots')
  dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
  output <- character()
  failures <- character()
  for (j in seq_len(nrow(targets))) {
    label <- as.character(targets$label[j])
    files <- list.files(run_dir,
      pattern=paste0('^', label, '_stratum_[0-9]+\\.rds$'), full.names=TRUE)
    if (!length(files)) {
      warning('No saved OSA residual objects for ', label)
      next
    }
    files <- sort(files)
    # Plot each stratum separately. Multiple run_osa objects may carry the same
    # fleet name, which can cause afscOSA::plot_osa() to fail when creating
    # factor levels for faceting. This changes only display labels, not residuals.
    for (k in seq_along(files)) {
      obj <- readRDS(files[k])
      stratum_id <- sprintf('stratum_%02d', k)
      display_name <- paste(label, stratum_id, sep='__')
      if (!is.null(obj$res) && 'fleet' %in% names(obj$res)) {
        obj$res$fleet <- rep(display_name, nrow(obj$res))
      }
      if (!is.null(obj$agg) && 'fleet' %in% names(obj$agg)) {
        obj$agg$fleet <- rep(display_name, nrow(obj$agg))
      }
      # Some package versions store additional tables carrying a fleet column.
      for (nm in names(obj)) {
        if (is.data.frame(obj[[nm]]) && 'fleet' %in% names(obj[[nm]])) {
          obj[[nm]]$fleet <- rep(display_name, nrow(obj[[nm]]))
        }
      }
      tryCatch({
        diagnostic <- afscOSA::plot_osa(list(obj),
                          use_agg_proportions=TRUE, plot=FALSE)
        dest <- file.path(plot_dir, paste0(display_name,'_afscOSA_composite.pdf'))
        .osa_composite_page(diagnostic, label=display_name,
                            destination=dest, show=show)
        output <- c(output, dest)
        saveRDS(diagnostic, file.path(plot_dir,
                         paste0(display_name,'_afscOSA_plots.rds')))
      }, error=function(e) {
        failures <<- c(failures, paste(display_name, conditionMessage(e)))
        warning('Could not plot ', display_name, ': ', conditionMessage(e),
                call.=FALSE)
      })
    }
  }
  if (length(failures)) {
    writeLines(failures, file.path(plot_dir,'plot_failures.txt'))
  }
  if (!length(output)) warning('No OSA diagnostic plots were created; inspect plot_failures.txt')
  output
}

osa_reweight_ss3 <- function(model_dir, ss_exe=NULL, max_iter=8L,
                             tolerance=0.05, damping=0.5, step_limit=1.6,
                             weight_bounds=c(0.02,3), seed=10L,
                             output_dir=NULL, min_years=3L,
                             dry_run=FALSE, plot_final=FALSE, show_final_plots=FALSE) {
  if(!requireNamespace('r4ss',quietly=TRUE) || !requireNamespace('afscOSA',quietly=TRUE))
    stop('Install packages r4ss and afscOSA first')
  if(!dir.exists(model_dir)) stop('Model folder does not exist: ',model_dir)
  model_dir <- normalizePath(model_dir)
  start <- .osa_starter(file.path(model_dir,'starter.ss'))
  ctl <- file.path(model_dir,start$control)
  dat <- file.path(model_dir,start$data)
  if(!file.exists(ctl) || !file.exists(dat)) stop('Missing data or control file named in starter.ss: ',start$data,' / ',start$control)
  if(max_iter<1 || damping<=0 || damping>1 || step_limit<1 || tolerance<=0 ||
     length(weight_bounds)!=2L || weight_bounds[1]<=0 || weight_bounds[2]<=weight_bounds[1]) stop('Invalid tuning arguments')
  rows <- .osa_rows(ctl)
  w <- rows[rows$factor %in% c(4L,5L) & rows$weight>0,,drop=FALSE]
  if(!nrow(w)) stop('No positive multinomial length (factor 4) or age (factor 5) variance adjustments found')
  w$type <- ifelse(w$factor==4L,'length','age')
  w$label <- paste0(w$type,'_fleet_',w$fleet)
  if(is.null(output_dir)) output_dir <- file.path(model_dir,'OSA_reweight_runs')
  if(dir.exists(output_dir)) stop('Output already exists; choose a new output_dir: ',output_dir)
  dir.create(output_dir,recursive=TRUE)
  output_dir <- normalizePath(output_dir)
  write.csv(w,file.path(output_dir,'detected_targets.csv'),row.names=FALSE)
  if(is.null(ss_exe)) {
    candidates <- c(file.path(model_dir,'ss3.exe'),file.path(model_dir,'ss3'),
                    file.path(model_dir,'SS3.exe'),Sys.which('ss3'),Sys.which('ss3.exe'))
    candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
    exe <- if(length(candidates)) normalizePath(candidates[1]) else ''
  } else {
    exe <- if(file.exists(ss_exe)) normalizePath(ss_exe) else Sys.which(ss_exe)
  }
  if(!nzchar(exe) && !dry_run) stop('SS3 executable not found. Supply ss_exe= with the executable path.')
  history <- list(); converged <- FALSE
  for(iter in seq_len(max_iter)) {
    run_dir <- file.path(output_dir,sprintf('iter_%02d',iter));dir.create(run_dir)
    .osa_copy_inputs(model_dir,run_dir,output_dir)
    .osa_write(ctl,file.path(run_dir,start$control),w)
    if(dry_run) { message('Prepared: ',run_dir); break }
    old <- getwd();setwd(run_dir)
    status <- tryCatch(system2(exe,args='-nohess',stdout='ss_stdout.txt',stderr='ss_stderr.txt',wait=TRUE),
                       finally=setwd(old))
    if(!is.numeric(status) || status!=0L || !file.exists(file.path(run_dir,'Report.sso')))
      stop('SS3 failed or did not write Report.sso in ',run_dir)
    fit <- r4ss::SS_output(dir=run_dir,verbose=FALSE,printstats=FALSE,covar=FALSE,
                           forecast=FALSE,warn=FALSE)
    records <- vector('list',nrow(w))
    for(j in seq_len(nrow(w))) {
      tbl <- if(w$factor[j]==4L) fit$lendbase else fit$agedbase
      if(is.null(tbl) || !nrow(tbl)) stop('Missing fitted composition table for ',w$label[j])
      parts <- .osa_prepare(as.data.frame(tbl),w$fleet[j],w$type[j],min_years=min_years)
      if(!length(parts)) stop('No usable composition strata for ',w$label[j],'; not safe to reweight')
      rr <- numeric(); metrics <- list()
      for(k in seq_along(parts)) {
        d <- parts[[k]]
        # Different random streams, reproducible across runs and iterations.
        call_seed <- as.integer(seed + 1000L*j + k)
        out <- afscOSA::run_osa(fleet=w$label[j],index_label=paste0('stratum_',k),
                 obs=d$obs,exp=d$exp,N=d$N,index=d$index,years=d$years,seed=call_seed)
        saveRDS(out,file.path(run_dir,sprintf('%s_stratum_%02d.rds',w$label[j],k)))
        z <- .osa_residuals(out)
        if(length(z)<2) stop('Insufficient OSA residuals in ',w$label[j],' stratum ',k)
        rr <- c(rr,z)
        metrics[[k]] <- data.frame(target=w$label[j],stratum=d$stratum,
                       n_years=length(d$years),n_resid=length(z),mean=mean(z),sdnr=sd(z))
      }
      write.csv(do.call(rbind,metrics),file.path(run_dir,paste0(w$label[j],'_strata.csv')),row.names=FALSE)
      if(length(rr)<10 || !is.finite(sd(rr)) || sd(rr)<=0) stop('Unusable pooled SDNR for ',w$label[j])
      records[[j]] <- data.frame(iter=iter,target=w$label[j],factor=w$factor[j],fleet=w$fleet[j],
                  old_weight=w$weight[j],SDNR=sd(rr),mean_resid=mean(rr),n_resid=length(rr))
    }
    result <- do.call(rbind,records)
    # Proposed damped heuristic: factor 1/SDNR^2, NOT established methodology.
    delta <- pmax(-log(step_limit),pmin(log(step_limit),-2*damping*log(result$SDNR)))
    result$proposed_next_weight <- pmax(weight_bounds[1],pmin(weight_bounds[2],result$old_weight*exp(delta)))
    result$within_tolerance <- abs(result$SDNR-1)<=tolerance
    write.csv(result,file.path(run_dir,'osa_weight_summary.csv'),row.names=FALSE)
    history[[iter]] <- result
    write.csv(do.call(rbind,history),file.path(output_dir,'osa_weight_history.csv'),row.names=FALSE)
    print(result[,c('target','old_weight','SDNR','proposed_next_weight','n_resid')],row.names=FALSE)
    if(all(result$within_tolerance)) { converged <- TRUE;break }
    if(iter<max_iter) w$weight <- result$proposed_next_weight
  }
  final_plot_files <- character()
  if (length(history) && isTRUE(plot_final)) {
    final_dir <- file.path(output_dir, sprintf('iter_%02d', length(history)))
    final_plot_files <- .osa_plot_final(final_dir, w, show=show_final_plots)
    message('Final afscOSA diagnostics: ', file.path(final_dir, 'final_osa_plots'))
  }
  message('Final fitted model: ',if(length(history)) file.path(output_dir,sprintf('iter_%02d',length(history))) else '(dry run)')
  message('Last proposed weights are NOT fitted unless they appear as old_weight in a subsequent iteration.')
  invisible(list(output_dir=output_dir,history=if(length(history)) do.call(rbind,history) else NULL,
                 converged=converged,targets=w,final_plot_files=final_plot_files))
}
