# TMLE-RMTL final analysis code
#
# This simplified version keeps the final estimator used in the manuscript:
# cause-specific Cox or HAL initial estimators, censoring and treatment weighting,
# one-step TMLE targeting of RMTL under competing risks, and EIC-based standard errors.
# Development-only switches and debugging parameters were removed from the public interface.

library(data.table)
library(survival)
library(Matrix)
library(zoo)
library(glmnet)
library(stringr)

lebesgue.loss.fun <- function(train.fit, dt, risk.set=NULL, test.set, time.var=NULL, X=NULL, lambda.cv=NULL,
                              delta.var="delta", delta.value=1) {

  tmp <- copy(dt)

  if (inherits(train.fit, "glmnet")) {
    beta.dim <- dim(train.fit$beta)
    if (length(beta.dim)==2 && (beta.dim[1]==0 || beta.dim[2]==0)) {
      return(Inf)
    }
  }

  pred.lambda <- try(predict(train.fit, X[tmp$id %in% test.set,],
                             newoffset=0, s=lambda.cv), silent=TRUE)
  if (inherits(pred.lambda, "try-error")) {
    return(Inf)
  }
  pred.lambda <- as.numeric(pred.lambda)
  if (length(pred.lambda) != sum(tmp$id %in% test.set) || any(!is.finite(pred.lambda))) {
    return(Inf)
  }

  tmp[id %in% test.set, fit.lambda:=exp(pred.lambda)]
  tmp[id %in% test.set, fit.dLambda:=fit.lambda*risk.time]
  tmp[id %in% test.set, fit.Lambda:=cumsum(fit.dLambda), by="id"]

  tmp[id %in% test.set, dN:=1*(time.obs==grid.time & get(delta.var)==delta.value)]

  return(tmp[id %in% test.set & time.obs==grid.time, -sum(log(fit.lambda)*dN - fit.Lambda)])
}
hal.screening <- function(covars, dt, V=5, cut.one.way=18, method.risk="test", cv.glmnet=FALSE,
                          seed=13349, grouped=TRUE, use.min=TRUE, order=1,
                          cut.time.treatment=NULL,
                          delta.var="delta", delta.value=1,
                          treatment=NULL, time.var=NULL, cut.time=5,
                          mat=NULL, browse=FALSE,
                          cut.two.way=5) {

  if (order==1) {
    fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way, method.risk=method.risk,
                   cut.time.treatment=cut.time.treatment,
                   delta.var=delta.var, delta.value=delta.value,
                   mat=mat, browse=browse,
                   cv.glmnet=cv.glmnet, treatment=treatment, time.var=time.var, cut.time=cut.time,
                   seed=seed, grouped=grouped, use.min=use.min)
    picked.covars <- covars[sapply(covars, function(covar) length(grep(covar, coef(fit$fit)@Dimnames[[1]][coef(fit$fit)@i+1]))>0)]
    if (length(picked.covars)==0) {
      for (jj in 1:10) {
        fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way*2, method.risk=method.risk,
                       cut.time.treatment=cut.time.treatment,
                       delta.var=delta.var, delta.value=delta.value,
                       mat=mat, browse=browse, lambda.cvs=fit$cve$min$lambda.cv/3,
                       cv.glmnet=cv.glmnet, treatment=treatment, time.var=time.var, cut.time=cut.time,
                       seed=seed, grouped=grouped, use.min=use.min)
        picked.covars <- covars[sapply(covars, function(covar) length(grep(covar, coef(fit$fit)@Dimnames[[1]][coef(fit$fit)@i+1]))>0)]
        if (length(picked.covars)>0) break
      }
    }
    return(picked.covars)
  } else {
    two.way <- expand.grid(c(treatment, covars), c(treatment, covars), stringsAsFactors=FALSE)
    fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way, method.risk=method.risk, cv.glmnet=cv.glmnet,
                   seed=seed, grouped=grouped, use.min=use.min, treatment=treatment,
                   cut.time.treatment=cut.time.treatment,
                   delta.var=delta.var, delta.value=delta.value,
                   mat=mat, browse=browse,
                   time.var=time.var, cut.time=cut.time,
                   two.way=two.way, cut.two.way=cut.two.way)$fit
    two.way.out <- two.way[apply(two.way, 1, function(row2) {
      if (row2[1]!=row2[2]) {
        (length(grep(row2[1], grep(paste0(":", row2[2]), coef(fit)@Dimnames[[1]][coef(fit)@i+1], value=TRUE)))>0)
      } else return(FALSE)
    }),, drop=FALSE]
    if (is.na(two.way.out[1,1])) two.way.out <- cbind(var1="", var2="")
    return(two.way.out)
  }

}
basis.fun <- function(covars, cut.one.way, dt, treatment=NULL,
                      Y=paste0("(", delta.var, "==", delta.value, ")"),
                      delta.var="delta", delta.value=1,
                      two.way=cbind(var1="", var2=""), cut.two.way=5,
                      cut.time.treatment=NULL, predict=FALSE, pseudo.dt=NULL,
                      time.var=NULL, cut.time=5) {

  if (length(covars)>0) {
    dt[, (covars):=lapply(.SD, function(x) {
      if (is.character(x)) {
        return(as.numeric(as.factor(x)))
      } else return(x)
    }), .SDcols=covars]
  }

  if (length(time.var)>0) {
    if (length(pseudo.dt)>0) {
      grid.times <- c(0, indicator.basis.fun(dt, time.var, cut.time, return.grid=TRUE), Inf)
      if (predict) grid.times <- sort(c(predict, grid.times))
      pseudo.dt[, grid.period:=as.numeric(cut(get(time.var), grid.times, include.lowest=TRUE, right=FALSE))]
      pseudo.dt <- unique(pseudo.dt, by=c("id", "grid.period"))
      pseudo.dt[, grid.time:=grid.times[grid.period]]
      if (!predict) pseudo.dt <- pseudo.dt[grid.time<time.obs]
      pseudo.dt[, grid.time:=grid.times[grid.period+1]]
      if (!predict) pseudo.dt[time.obs<=grid.time, grid.time:=time.obs]
      pseudo.dt[, (time.var):=grid.time]
    } else {
      pseudo.dt <- copy(dt)
      if (length(time.var)>0) {
        grid.times <- c(0, indicator.basis.fun(dt, time.var, cut.time, return.grid=TRUE), Inf)
        if (predict) grid.times <- sort(c(predict, grid.times))
        pseudo.dt <- do.call("rbind", lapply(1:length(grid.times), function(jj) {
          tmp <- copy(pseudo.dt)[, grid.period:=jj]
        }))
      }
      pseudo.dt[, grid.time:=grid.times[grid.period]]
      if (!predict) pseudo.dt <- pseudo.dt[grid.time<get(time.var)]
      pseudo.dt <- pseudo.dt[order(id, get(time.var))]
      pseudo.dt[, grid.time:=grid.times[grid.period+1]]
      if (!predict) pseudo.dt[get(time.var)<=grid.time, grid.time:=get(time.var)]
      pseudo.dt[, time.obs:=get(time.var)]
      pseudo.dt[, (time.var):=grid.time]
    }
    pseudo.dt <- pseudo.dt[!is.na(time)]
  } else {
    pseudo.dt <- copy(dt)
  }

  X <- Matrix(
    model.matrix(formula(paste0(
      Y, "~-1",
      ifelse(length(treatment)>0, paste0("+", treatment), ""),
      ifelse(length(time.var)>0, paste0("+",
                                        paste0(indicator.basis.fun(pseudo.dt, "grid.period", cut.time), collapse="+")), ""),
      ifelse(length(time.var)>0 & cut.time.treatment>0,
             paste0("+", paste0(paste0(indicator.basis.fun(pseudo.dt, "grid.period", cut.time.treatment), ":", treatment, "+"), collapse="")), ""),
      ifelse(length(covars)>0 & cut.one.way>2,
             paste0("+", paste0(sapply(covars, function(covar) paste0(indicator.basis.fun(pseudo.dt, covar, cut.one.way), collapse="+")), collapse="+")),
             ""),
      ifelse(two.way[1,1]!="" & cut.two.way>2,
             paste0(apply(two.way, 1, function(row) {
               if (row[1]!=row[2]) {
                 if (any(row==treatment)) {
                   paste0("+", paste0(apply(expand.grid(ifelse(row[1]==treatment, treatment, indicator.basis.fun(dt, row[1], cut.two.way)),
                                                        ifelse(row[2]==treatment, treatment, indicator.basis.fun(dt, row[2], cut.two.way))), 1,
                                            function(row2) paste0(row2, collapse=":")), collapse="+"))
                 } else {
                   paste0("+", paste0(apply(expand.grid(indicator.basis.fun(dt, row[1], cut.two.way),
                                                        indicator.basis.fun(dt, row[2], cut.two.way)), 1,
                                            function(row2) paste0(row2, collapse=":")), collapse="+"))
                 }
               } else return("")
             }), collapse=""), "")
    )),
    data=pseudo.dt), sparse=FALSE)

  if (length(time.var)>0) {
    x.vector <- apply(X, 1, function(x) paste0(x, collapse=","))
    pseudo.dt[, x:=x.vector]
    pseudo.dt[, risk.time:=grid.time-c(0, grid.time[-.N]), by="id"]
    return(list(X=X, pseudo.dt=pseudo.dt))
  } else {
    return(X)
  }

}

indicator.basis.fun <- function(mat, xvar, xcut, type="obs", return.grid=FALSE, seed=220) {
  set.seed(seed)
  if (type=="obs") {
    xvar.values <- mat[, sort(unique(get(xvar)))]
    xvar.pick <- seq(1, length(xvar.values), length=min(xcut, length(xvar.values)))
    xgrid <- xvar.values[xvar.pick][-c(1,xcut)]
  } else {
    xgrid <- round(seq(mat[, min(get(xvar))],
                       mat[, max(get(xvar))],
                       length=xcut)[-c(1,xcut)], 2)
  }
    if (return.grid) return(c(unique(xgrid))) else return(paste0("(", xvar, ">=", unique(xgrid), ")"))
}

align.hal.design <- function(X.new, X.ref) {
    ref.names <- colnames(X.ref)
    new.names <- colnames(X.new)

    X.new <- as.matrix(X.new)
    X.ref <- as.matrix(X.ref)

    if (is.null(ref.names) | is.null(new.names)) {
        if (ncol(X.new) == ncol(X.ref)) return(X.new)
        stop("HAL design matrices could not be aligned because column names were missing.")
    }

    missing.names <- setdiff(ref.names, new.names)
    if (length(missing.names) > 0) {
        add.mat <- matrix(0, nrow=nrow(X.new), ncol=length(missing.names))
        colnames(add.mat) <- missing.names
        X.new <- cbind(X.new, add.mat)
    }

    extra.names <- setdiff(colnames(X.new), ref.names)
    if (length(extra.names) > 0) {
        X.new <- X.new[, !(colnames(X.new) %in% extra.names), drop=FALSE]
    }

    X.new[, ref.names, drop=FALSE]
}
cv.fun <- function(loss.fun, dt, V=5, seed=19192, X=NULL, Y=NULL, offset=NULL,
                   method.risk=c("test","train","VvH"), time.var=NULL,
                   penalty.factor=rep(1, ncol(X)), delta.var="delta", delta.value=1,
                   change.point=NULL, treatment=NULL,
                   verbose=FALSE,
                   cox.model=NULL, lambda.cvs=c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))) {

  set.seed(seed)

  n <- length(dt[, unique(id)])

  cv <- sample(1:n, size=n)
  V2 <- ceiling(n/V)
  cv.split <- lapply(1:V, function(vv) na.omit(cv[((vv-1)*V2+1):(vv*V2+1)]))

  cve.list <- lapply(1:V, function(vv) {

    test.set <- cv.split[[vv]]
    train.set <- dt[, id][!dt[, id] %in% test.set]

    dt.train <- dt[id %in% train.set]

    if (length(X)==0) {

      if (ifelse(length(change.point)>0, change.point>0, FALSE)) {

        dt2.train <- rbind(dt.train, dt.train)[order(id)]
        dt2.train[, time.period:=1:.N, by="id"]
        dt2.train[, time.indicator:=(time<=change.point)]
        dt2.train[time.period==1, `:=`(tstart=0, tstop=(time<=change.point)*time+(time>change.point)*change.point)]
        dt2.train[time.period==2, `:=`(tstart=change.point, tstop=time)]
        dt2.train[time.period==1 & !time.indicator, (delta.var):=delta.value-1]
        dt2.train <- dt2.train[!(time.period==2 & time.indicator)]

        mod1 <- as.character(cox.model)
        mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                   which(strsplit(mod1[2], "")[[1]]==",")-1), paste0("tstart", ", tstop"),
                            mod1[2]),
                       "~",
                       gsub("\\+\\+", "", gsub(paste0("\\+", treatment), "",
                                               gsub(paste0("\\+", treatment), "",
                                                    gsub(" ", "",
                                                         paste0("I((time.period","==1)&(", treatment, "==1))",
                                                                " + I((time.period", "==2)&(", treatment, "==1))",
                                                                paste0("+", mod1[3])))))))
        train.fit <- coxph(as.formula(mod2), data=dt2.train)

      } else {
        train.fit <- coxph(cox.model, data=dt.train)
      }

      if (str_to_lower(method.risk[1])=="vvh") {
        return(loss.fun(train.fit=train.fit, dt=dt, risk.set=1:n, test.set=test.set, delta.var=delta.var, delta.value=delta.value, change.point=change.point)-
                 loss.fun(train.fit=train.fit, dt=dt, risk.set=train.set, test.set=test.set, delta.var=delta.var, delta.value=delta.value, change.point=change.point))
      } else if (method.risk[1]=="test") {
        return(loss.fun(train.fit=train.fit, dt=dt, risk.set=test.set, test.set=test.set, delta.var=delta.var, delta.value=delta.value, change.point=change.point))
      } else {
        return(loss.fun(train.fit=train.fit, dt=dt, risk.set=train.set, test.set=test.set, delta.var=delta.var, delta.value=delta.value, change.point=change.point))
      }

    } else {

      if (length(offset)>0) {

        dt.train <- dt[id %in% train.set]

        dt.train[, D:=sum(time.obs==grid.time & get(delta.var)==delta.value), by="x"]
        dt.train[, risk.time:=grid.time-c(0, grid.time[-.N]), by="id"]
        dt.train[, RT:=sum(risk.time), by="x"]

        tmp.dt <- unique(dt.train[, c("x", "D", "RT")])
        Y.train <- tmp.dt[RT>0, D]
        offset.train <- tmp.dt[RT>0, log(RT)]
        X.train <- unique.matrix(X[dt$id %in% train.set,])[tmp.dt$RT>0,]

        if (length(lambda.cvs)>0) {

          train.fit <- glmnet(x=as.matrix(X.train), y=Y.train,
                              offset=offset.train,
                              family="poisson",
                              maxit=1000,
                              penalty.factor=penalty.factor,
                              lambda=lambda.cvs)

        } else {

          train.fit <- cv.glmnet(x=as.matrix(X.train), y=Y.train,
                                 offset=offset.train,
                                 family="poisson",
                                 penalty.factor=penalty.factor,
                                 maxit=1000)

          lambda.cvs <- train.fit$lambda.1se

        }

        return(sapply(lambda.cvs, function(lambda.cv) {
          return(loss.fun(train.fit=train.fit, dt=dt, test.set=test.set,
                          time.var=time.var,
                          X=X, lambda.cv=lambda.cv, delta.var=delta.var, delta.value=delta.value))
        }))

      } else {

        Y.train <- Y[dt$id %in% train.set]
        X.train <- X[dt$id %in% train.set,]

        if (length(lambda.cvs)>0) {

          train.fit <- glmnet(x=as.matrix(X.train), y=Y.train,
                              family="cox",
                              maxit=1000,
                              penalty.factor=penalty.factor,
                              lambda=lambda.cvs)

        } else {

          train.fit <- cv.glmnet(x=as.matrix(X.train), y=Y.train,
                                 family="cox",
                                 penalty.factor=penalty.factor,
                                 maxit=1000)

          lambda.cvs <- train.fit$lambda.1se

        }

        return(sapply(lambda.cvs, function(lambda.cv)
          if (str_to_lower(method.risk[1])=="vvh") {
            loss.fun(train.fit=train.fit, dt=dt, risk.set=1:n, test.set=test.set,
                     X=X, lambda.cv=lambda.cv, delta.var=delta.var, delta.value=delta.value)-
              loss.fun(train.fit=train.fit, dt=dt, risk.set=train.set, test.set=test.set,
                       X=X, lambda.cv=lambda.cv, delta.var=delta.var, delta.value=delta.value)
          } else if (method.risk[1]=="test") {
            loss.fun(train.fit=train.fit, dt=dt, risk.set=test.set, test.set=test.set,
                     X=X, lambda.cv=lambda.cv, delta.var=delta.var, delta.value=delta.value)
          } else {
            loss.fun(train.fit=train.fit, dt=dt, risk.set=train.set, test.set=test.set,
                     X=X, lambda.cv=lambda.cv, delta.var=delta.var, delta.value=delta.value)
          }))
      }

    }
  })

  if (length(X)==0 | length(lambda.cvs)==0) {
    return(list(cve=sum(unlist(cve.list))))
  } else {
    cve <- unlist(lapply(1:length(lambda.cvs), function(mm) {
      sum(unlist(lapply(cve.list, function(out) out[[mm]])))
    }))
    lambda.cv <- min(lambda.cvs[cve==min(cve)])
    return(list(min=list(lambda.cv=lambda.cv,
                         cve=unique(cve[cve==min(cve)])),
                all=cbind(lambda=lambda.cvs, cve=cve)))
  }

}
cox.loss.fun <- function(train.fit, dt, risk.set, test.set, X=NULL, lambda.cv=NULL,
                         delta.var="delta", delta.value=1, change.point=NULL) {

  tmp <- copy(dt)

  if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
    tmp <- rbind(tmp, tmp)[order(id)]
    tmp[, time.period:=1:.N, by="id"]
  }

  if (any(class(train.fit)%in%c("coxnet", "cv.glmnet", "glmnet"))) {
    tmp[, fit.lp:=predict(train.fit, type="link", newx=X, s=lambda.cv)]
  } else {
    tmp[, fit.lp:=predict(train.fit, type="lp", newdata=tmp)]
  }

  tmp[, risk:=0]
  tmp[id%in%risk.set, risk:=1]

  tmp <- tmp[rev(order(time))]

  if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
    tmp[, term2:=cumsum(risk*exp(fit.lp)), by="time.period"]
  } else {
    tmp[, term2:=cumsum(risk*exp(fit.lp))]
  }

  tmp[term2==0, term2:=1]

  if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
    return(-sum(tmp[id%in%test.set &
                      ((time.period==1 & time<=change.point) | (time.period==2 & time>change.point)),
                    (get(delta.var)==delta.value)*(fit.lp - log(term2))]))
  } else {
    return(-sum(tmp[id%in%test.set, (get(delta.var)==delta.value)*(fit.lp - log(term2))]))
  }

}
fit.hal <- function(covars, dt, treatment=NULL, V=5, cut.one.way=8, method.risk="test", cv.glmnet=FALSE,
                    seed=13349, grouped=TRUE, use.min=TRUE, penalize.treatment=FALSE,
                    penalize.time=TRUE, treatment.prediction=treatment,
                    intervention=c(0,1),
                    delta.var="delta", delta.value=1,
                    time.var=NULL, cut.time=5, cut.time.treatment=NULL,
                    browse=FALSE, predict=FALSE, mat=NULL, verbose=FALSE,
                    lambda.cvs=c(sapply(1:5, function(jjj) (9:1)/(10^jjj))),
                    two.way=cbind(var1="", var2=""), cut.two.way=5,
                    return.cve=FALSE) {

  set.seed(seed)

  X.hal <- basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=dt,
                     time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                     delta.var=delta.var, delta.value=delta.value, pseudo.dt=mat,
                     two.way=two.way, cut.two.way=cut.two.way)

  if (is.list(X.hal) & length(X.hal)==2) {
    pseudo.dt <- X.hal$pseudo.dt
    X.hal <- X.hal$X
    pseudo.dt[, D:=sum(time.obs==grid.time & get(delta.var)==delta.value), by="x"]
    pseudo.dt[, RT:=sum(risk.time), by="x"]
  } else {
    Y.hal <- dt[, Surv(time, get(delta.var)==delta.value)]
  }

  penalty.factor <- rep(1, ncol(X.hal))

  if (!penalize.treatment & length(grep(treatment, colnames(X.hal)))>0) {
    penalty.factor[colnames(X.hal) %in% grep(treatment, colnames(X.hal), value=TRUE)] <- 0
    penalty.factor <- 1-(cumprod(1-penalty.factor))
  }

  if (!penalize.time & length(grep("grid.period", colnames(X.hal)))>0) {
    not.penalize <- grep("grid.period", colnames(X.hal), value=TRUE)
    if (length(grep(":", not.penalize))>0) not.penalize <- not.penalize[-grep(":", not.penalize)]
    penalty.factor[colnames(X.hal) %in% not.penalize] <- 0
  }

  if (length(time.var)>0) {
    tmp.dt <- unique(pseudo.dt[, c("x", "D", "RT")])
    Y2.hal <- tmp.dt[RT>0, D]
    offset2 <- tmp.dt[RT>0, log(RT)]
    X2.hal <- unique.matrix(X.hal)[tmp.dt$RT>0,]
    if (cv.glmnet) {
      hal.fit <- cv.glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                           offset=offset2,
                           family="poisson",
                           penalty.factor=penalty.factor,
                           maxit=1000)
      if (use.min) {
        lambda.cv <- hal.fit$lambda.min
        hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                          offset=offset2,
                          family="poisson",
                          lambda=lambda.cv,
                          penalty.factor=penalty.factor,
                          maxit=1000)
        cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                          penalty.factor=penalty.factor,
                          offset=TRUE, V=V,lambda.cvs=lambda.cv, delta.var=delta.var, delta.value=delta.value)
      } else {
        cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                          penalty.factor=penalty.factor,
                          offset=TRUE, V=V,lambda.cvs=hal.fit$lambda.1se, delta.var=delta.var, delta.value=delta.value)
      }
    } else {
      cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                        penalty.factor=penalty.factor, lambda.cvs=lambda.cvs,
                        offset=TRUE, V=V, delta.var=delta.var, delta.value=delta.value)
      hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                        offset=offset2,,
                        lambda=cve.hal$min$lambda.cv,
                        family="poisson",
                        penalty.factor=penalty.factor,
                        maxit=1000)
      repeat {
        if (sum(coef(hal.fit))==0) {
          cve.hal <- cv.fun(loss.fun=lebesgue.loss.fun, dt=pseudo.dt, X=X.hal, Y=NULL, time.var=time.var,
                            penalty.factor=penalty.factor, lambda.cvs=lambda.cvs[lambda.cvs>cve.hal$min$lambda.cv],
                            offset=TRUE, V=V, delta.var=delta.var, delta.value=delta.value)
          hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                            offset=offset2,,
                            lambda=cve.hal$min$lambda.cv,
                            family="poisson",
                            penalty.factor=penalty.factor,
                            maxit=1000)
        } else {
          break
        }
      }
      lambda.cv <- cve.hal$min$lambda.cv
      if (browse) {
        lambda.cvs <- c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))
        hal.fit <- glmnet(x=as.matrix(X2.hal), y=Y2.hal,
                          offset=offset2,,
                          lambda=lambda.cvs,
                          family="poisson",
                          penalty.factor=penalty.factor,
                          maxit=1000)
        lapply(lambda.cvs, function(lambda.cv) covars[sapply(covars, function(covar) length(grep(covar, coef(hal.fit, s=lambda.cv)@Dimnames[[1]][coef(hal.fit, s=lambda.cv)@i+1]))>0)])
      }
    }
  } else if (cv.glmnet) {
    hal.fit <- cv.glmnet(x=as.matrix(X.hal), y=Y.hal,
                         family="cox",
                         grouped=grouped,
                         penalty.factor=penalty.factor,
                         maxit=1000)
    if (use.min) {
      lambda.cv <- hal.fit$lambda.min
      hal.fit <- glmnet(x=as.matrix(X.hal), y=Y.hal,
                        family="cox",
                        lambda=lambda.cv,
                        grouped=grouped,
                        penalty.factor=penalty.factor,
                        maxit=1000)
      cve.hal <- cv.fun(loss.fun=cox.loss.fun, dt=dt, X=X.hal, Y=Y.hal, method.risk=method.risk, V=V,
                        penalty.factor=penalty.factor,
                        lambda.cvs=lambda.cv, delta.var=delta.var, delta.value=delta.value)
    } else {
      cve.hal <- cv.fun(loss.fun=cox.loss.fun, dt=dt, X=X.hal, Y=Y.hal, method.risk=method.risk, V=V,
                        penalty.factor=penalty.factor,
                        lambda.cvs=hal.fit$lambda.1se, delta.var=delta.var, delta.value=delta.value)
    }
  } else {
    cve.hal <- cv.fun(loss.fun=cox.loss.fun, dt=dt, X=X.hal, Y=Y.hal, lambda.cvs=lambda.cvs,
                      method.risk=method.risk, V=V, delta.var=delta.var, delta.value=delta.value)
    hal.fit <- glmnet(x=as.matrix(X.hal), y=Y.hal,
                      lambda=cve.hal$min$lambda.cv,
                      family="cox",
                      penalty.factor=penalty.factor,
                      maxit=1000)
    lambda.cv <- cve.hal$min$lambda.cv
    if (browse) {
      warning("Interactive HAL browsing is disabled in the simplified final code.")
      lambda.cvs <- c(sapply(1:5, function(jjj) (9:1)/(10^jjj)))
      hal.fit <- glmnet(x=as.matrix(X.hal), y=Y.hal,
                        lambda=cve.hal$min$lambda.cvs,
                        family="cox",
                        penalty.factor=penalty.factor,
                        maxit=1000)
      lapply(lambda.cvs, function(lambda.cv) covars[sapply(covars, function(covar) length(grep(covar, coef(hal.fit, s=lambda.cv)@Dimnames[[1]][coef(hal.fit, s=lambda.cv)@i+1]))>0)])
    }
  }

  if (verbose) {
    print(coef(hal.fit, s=lambda.cv))
  }

  if (predict) {

    if (length(time.var)>0) {

      if (length(mat)>0) {
        X.hal.a <- lapply(intervention, function(a) {
          return(basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=a],
                           time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                           predict=predict,
                           pseudo.dt=copy(mat)[,(treatment):=a],
                           delta.var=delta.var, delta.value=delta.value,
                           two.way=two.way, cut.two.way=cut.two.way))
        })
      } else {
        X.hal.a <- lapply(intervention, function(a) {
          return(basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way, dt=copy(dt)[,(treatment):=a],
                           time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                           predict=predict,
                           delta.var=delta.var, delta.value=delta.value,
                           two.way=two.way, cut.two.way=cut.two.way))
        })
      }

      pseudo.dt <- X.hal.a[[1]]$pseudo.dt
      X.hal.a <- lapply(1:length(intervention), function(aa) {
        align.hal.design(X.hal.a[[aa]]$X, X.hal)
      })

      for (aa in 1:length(intervention)) {
        pseudo.dt[, (paste0("fit.lambda", intervention[aa])):=exp(predict(hal.fit, X.hal.a[[aa]],
                                                                          newoffset=0, s=lambda.cv))]
      }

      if (length(mat)>0) {

        mat <- do.call("rbind", lapply(intervention, function(a) {
          tmp_dt <- merge(mat[get(treatment.prediction)==a],
                          pseudo.dt[, c("id", covars, "grid.period", paste0("fit.lambda", a)), with=FALSE],
                          by=c("id", covars, "grid.period"),
                          all.x=TRUE)
          tmp_dt[, fit.lambda:=get(paste0("fit.lambda", a))]
          tmp_dt[, (paste0("fit.lambda", a)):=NULL]
          return(tmp_dt)
        }))

        setorderv(mat, c("id", treatment.prediction, time.var))
        mat[, risk.time:=c(0, diff(get(time.var))), by=c("id", treatment.prediction)]
        if (mat[is.na(fit.lambda) & risk.time>0, .N] > 0) {
          warning("Poisson-HAL prediction was missing for positive-risk intervals; setting those increments to zero.")
        }
        mat[is.na(fit.lambda), fit.lambda:=0]
        mat[, fit.dLambda:=fit.lambda*risk.time]
        mat[, fit.Lambda:=cumsum(fit.dLambda), by=c("id", treatment.prediction)]

        if (delta.value>0) {
          mat[, (paste0("dhaz", delta.value)):=fit.dLambda]
          mat[, (paste0("fit.cox", delta.value)):=1]
        } else {
          mat[, surv.C.pois:=exp(-cumsum(fit.dLambda)), by=c("id", treatment.prediction)]
          mat[, surv.C1.pois:=c(1, surv.C.pois[-.N]), by=c("id", treatment.prediction)]
          mat[, Ht:=Ht*surv.C1/surv.C1.pois]
          mat[, surv.C1:=surv.C1.pois]
        }

        if (return.cve) {
          return(list(mat=mat, cve=cve.hal))
        } else {
          return(mat)
        }

      } else {

        for (a in intervention) {
          pseudo.dt[, (paste0("fit.dLambda", a)):=get(paste0("fit.lambda", a))*risk.time]
          pseudo.dt[, (paste0("fit.Lambda", a)):=cumsum(get((paste0("fit.dLambda", a)))), by="id"]
          pseudo.dt[, (paste0("surv", a)):=exp(-get(paste0("fit.Lambda", a))), by="id"]
        }

        if (length(intervention)>1) {
          return(list(fit=hal.fit,
                      cve=cve.hal,
                      target.fit=pseudo.dt[time==predict, -mean(surv1 - surv0)]))
        } else {
          return(list(fit=hal.fit,
                      cve=cve.hal,
                      target.fit=pseudo.dt[time==predict, -1-mean(get(paste0("surv", intervention)))]))
        }
      }

    } else {

      X.hal.a <- lapply(intervention, function(a) {
        return(basis.fun(covars=covars, treatment=treatment, cut.one.way=cut.one.way,
                         dt=copy(dt)[,(treatment):=a],
                         time.var=time.var, cut.time=cut.time, cut.time.treatment=cut.time.treatment,
                         predict=predict,
                         delta.var=delta.var, delta.value=delta.value,
                         two.way=two.way, cut.two.way=cut.two.way))
      })
      X.hal.a <- lapply(X.hal.a, align.hal.design, X.ref=X.hal)

      survfit.hal <- lapply(1:length(intervention), function(aa) {
        survival::survfit(hal.fit, x=X.hal, y=Y.hal, s=lambda.cv,
                          newx=X.hal.a[[aa]])})

      predict.hal <- do.call("rbind", lapply(1:length(intervention), function(aa) {
        return(do.call("rbind", lapply(1:ncol(survfit.hal[[aa]]$surv), function(jj) {
          tmp <- data.table(id=jj,
                            time=survfit.hal[[aa]]$time,
                            surv=survfit.hal[[aa]]$surv[,jj],
                            cumhaz=survfit.hal[[aa]]$cumhaz[,jj])
          tmp[, (treatment):=intervention[aa]]
          if (!0 %in% tmp$time) {
            row0 <- data.table(id=jj, time=0, surv=1, cumhaz=0)
            row0[, (treatment):=intervention[aa]]
            tmp <- rbind(row0, tmp)
          }
          setorder(tmp, time)
          tmp[, dhaz:=c(cumhaz[1], diff(cumhaz))]
          return(tmp)
        })))
      }))

      if (length(mat)>0) {

        if ("cumhaz" %in% names(mat))
          mat[, cumhaz:=NULL]

        if ("dhaz" %in% names(mat))
          mat[, dhaz:=NULL]

        mat <- merge(mat, predict.hal, by=c("id", "time", treatment), all.x = TRUE)

        time.order.var <- if (length(time.var)>0) time.var else "time"
        setorderv(mat, c("id", treatment, time.order.var))
        if ("dhaz" %in% names(mat)) mat[is.na(dhaz), dhaz:=0]
        if ("cumhaz" %in% names(mat)) mat[, cumhaz:=cumsum(dhaz), by=c("id", treatment)]
        if ("surv" %in% names(mat)) mat[, surv:=exp(-cumhaz)]

        mat[, (paste0("fit.cox", delta.value)):=1]
        mat[, (paste0("chaz", delta.value)):=cumhaz]
        mat[, (paste0("dhaz", delta.value)):=dhaz]

        if (delta.value==0) {
          mat[, "surv.C.hal":=
                exp(-get(paste0("chaz", delta.value))*
                      get(paste0("fit.cox", delta.value)))]
          mat[, "surv.C1.hal":=c(1, surv.C.hal[-.N]), by=c("id", treatment)]
          mat[, Ht:=Ht*surv.C1/surv.C1.hal]
          mat[, surv.C1:=surv.C1.hal]
        }

        if (return.cve) {
          return(list(mat=mat, cve=cve.hal))
        } else {
          return(mat)
        }

      }

      if (length(intervention)>1) {
        warning("option to predict directly not maintained")
        return(list(fit=hal.fit,
                    cve=cve.hal,
                    target.fit=-mean(exp(-exp(predict(hal.fit, type="link", newx=X.hal1, s=lambda.cv))*basehaz[nrow(basehaz), hazard1])-
                                       exp(-exp(predict(hal.fit, type="link", newx=X.hal0, s=lambda.cv))*basehaz[nrow(basehaz), hazard0]))))
      } else {
        return(list(fit=hal.fit,
                    cve=cve.hal,
                    target.fit=1-mean(exp(-exp(predict(hal.fit, type="link", newx=X.hal.a[[1]], s=lambda.cv))*basehaz[nrow(basehaz), get(paste0("hazard", intervention))]))))
      }

    }

  }

  return(list(fit=hal.fit,
              cve=cve.hal))
}
contmle.rmtl <- function(
  dt,
  estimation = list(
    event1 = list(
      fit = "cox",
      model = Surv(time, delta == 1) ~ A + L1 + L2 + L3,
      changepoint = NULL
    ),
    event2 = list(
      fit = "cox",
      model = Surv(time, delta == 2) ~ A + L1 + L2 + L3,
      changepoint = NULL
    ),
    cens = list(
      fit = "cox",
      model = Surv(time, delta == 0) ~ A + L1 + L2 + L3,
      changepoint = NULL
    )
  ),
  target = 1,
  do.hal.screening = FALSE,
  hal.tuning = FALSE,
  treat.model = A ~ L1 + L2 + L3,
  treat.effect = c("ate", "1", "0"),
  tau = c(1.2),
  cut.off.cens = 1e-2,
  deps.size = 0.005,
  no.small.steps = 200,
  V = 10,
  lambda.cvs = seq(0.0000001, 0.01, length = 50),
  cut.one.way = 8,
  cut.time = 10,
  cut.time.A = 10,
  cut.two.way = 5,
  verbose = FALSE
) {
  hal.sl <- isTRUE(hal.tuning)
  sl.models <- list()
  verbose.hal <- verbose
  verbose.sl <- verbose
  weighted.norm <- FALSE
  eic.scale <- "current"
  trace.eic.scales <- FALSE
  no.update.if.solved <- FALSE
  simultaneous.ci <- FALSE
  use.observed.times <- FALSE
  length.times <- length(tau)
  output.mat <- FALSE
  cut.one.way.final <- cut.one.way
  cut.one.way.screen <- min(5, cut.one.way.final)
  cut.two.way.screen <- cut.two.way

  dt[, id:=1:.N]
  target.S=FALSE
  cr=TRUE
  weighted.norm <- weighted.norm[1]
  eic.scale <- eic.scale[1]
  if (!eic.scale %in% c("current", "initial")) {
    stop("eic.scale must be one of: current, initial")
  }
  if (is.null(cut.one.way.final)) cut.one.way.final <- cut.one.way
  if (is.null(cut.one.way.screen)) cut.one.way.screen <- min(5, cut.one.way.final)
  if (is.null(cut.two.way.screen)) cut.two.way.screen <- cut.two.way
  if (verbose) {
    verbose.sl <- TRUE
    verbose.hal <- TRUE
  }

  not.fit.list <- list()

  time.var <- gsub("Surv\\(", "", unlist(strsplit(as.character(estimation[[1]][["model"]])[2], ","))[1])
  delta.var <- gsub(" ", "",
                    gsub("==", "",
                         gsub("[0-9]+\\)", "",
                              unlist(strsplit(as.character(estimation[[1]][["model"]])[2], ","))[2])))

  estimation <- lapply(estimation, function(x) {
    event <- as.numeric(gsub("\\D", "", unlist(strsplit(as.character(x[["model"]])[2], ","))[2]))
    x[[length(x)+1]] <- event
    names(x)[length(x)] <- "event"
    if (!(event %in% dt[, unique(get(delta.var))]))
      message(paste0("model specified for ", delta.var, "=", event,
                     ", but there were no observations")) else return(x)
  })

  estimation <- estimation[sapply(estimation, function(each) length(each)>0)]

  outcome.index <- unlist(lapply(1:length(estimation), function(each) {
    event <- estimation[[each]][["event"]]
    each[event>0]
  }))

  events <- unlist(lapply(estimation, function(each) each[["event"]]))
  if (any(table(events)>1)) {
    stop(paste0("multiple models specified for ", delta.var, "=",
                paste0(names(table(events))[table(events)>1], collapse=",")))
  }

  delta.missing <- dt[, unique(get(delta.var))][!(dt[, unique(get(delta.var))] %in% unlist(lapply(estimation, function(each) each[["event"]])))]
  if (length(delta.missing)) {
    warning(paste0("No model specified for event type=", paste0(delta.missing, collapse=","),
                   "; will use the one specified for event type=",
                   estimation[[1]][["event"]]))
    for (delta in delta.missing) {
      estimation[[length(estimation)+1]] <- estimation[[1]]
      estimation[[length(estimation)]][["event"]] <- delta
      estimation[[length(estimation)]][["model"]] <-
        as.formula(paste0(gsub(estimation[[1]][["event"]], delta, estimation[[length(estimation)]][["model"]][2]),
                          estimation[[length(estimation)]][["model"]][1],
                          estimation[[length(estimation)]][["model"]][3]))
    }
  }

  n <- length(dt[, unique(id)])

  A.name <- as.character(treat.model)[2]

  covars <- NULL
  all_formulas <- c(lapply(sl.models, function(x) x[[1]]),
                    unlist(lapply(estimation, function(x) x[["model"]])))

  for (mod in all_formulas) {
    rhs <- as.character(mod)[3]
    vars <- unlist(strsplit(gsub("\\+", " ", rhs), " "))
    covars <- c(covars, vars)
  }

  covars <- unique(covars)
  covars <- covars[!covars %in% c(A.name, "", "*", "1",NA)]

  unique.times <- sort(unique(dt[, get(time.var)]))
  unique.times2 <- sort(unique(dt[get(delta.var)>0, get(time.var)]))

  if (use.observed.times) {
    n <- nrow(dt)
    unique.times3 <- unique.times[unique.times<=max(tau)]
    tau <- unique.times3[(1:length(unique.times3)) %% floor(length(unique.times3)/length.times)==1]
    sort(unique.times3[sample(length(unique.times3), length.times)])
  } else {
    test.tau <- findInterval(tau, unique.times2)
    if (FALSE & !length(test.tau)==length(unique(test.tau))) {
      tau <- na.omit(tau[(1:length(tau))[unique(findInterval(unique.times2, tau))+1]])
      warning("no observations between tau as specified, truncating tau")
    }
  }

  if (dt[, length(unique(get(A.name)))]>2) {
    a <- treat.effect
  } else {
    if (treat.effect[1]=="1") a <- 1 else if (treat.effect[1]=="0") a <- 0 else a <- c(1, 0)
  }

  dt2 <- NULL
  bhaz.cox <- do.call("rbind", lapply(a, function(aa) data.table(time=c(0, unique.times), A=aa)))
  setnames(bhaz.cox, "A", A.name)

  prob.A <- predict(glm(as.formula(treat.model), family = "binomial", data=dt), type="response")

  if (verbose & dt[, length(unique(get(A.name)))]==2) print(summary(glm(as.formula(treat.model),family = "binomial", data=dt)))
  if (verbose & dt[, length(unique(get(A.name)))]==2) print(summary(prob.A))

  for (each in 1:length(estimation)) {

    fit <- estimation[[each]][["fit"]][1]
    fit.model <- estimation[[each]][["model"]]
    if (any(names(estimation[[each]])=="changepoint"))
      fit.changepoint <- estimation[[each]][["changepoint"]] else fit.changepoint <- NULL
    fit.delta <- estimation[[each]][["event"]]
    fit.name <- names(estimation)[each]

    if (fit[1] %in% c("sl", "cox.hal.sl")) {
      stop("The simplified final code removes the development-only super-learner branch; use fit = 'cox', 'cox.hal', 'hal', or 'km'.")
    } else {

      sl.pick <- ""
      cve.sl.pick <- ""

      if (fit[1] %in% c("km")) {
        tmp.model <- as.character(fit.model)
        if (fit[1]=="km") tmp.model[3] <- "strata(A)" else tmp.model[3] <- "1"
        estimation[[each]]$model <- fit.model <- formula(paste0(tmp.model[2], tmp.model[1], tmp.model[3]))
        estimation[[each]]$changepoint <- fit.changepoint <- NULL
      }
    }

    estimation[[each]]$sl.pick <- fit.model
    estimation[[each]]$cve.sl.pick <- cve.sl.pick

    if (length(fit.changepoint)>0 & length(dt2)==0) {
      dt2 <- rbind(dt, dt)[order(id)]
    }

    if (length(fit.changepoint)>0) {
      delta1 <- abs(fit.delta-1)
      dt2[, time.indicator:=(get(time.var)<=fit.changepoint)]
      dt2[, (paste0("period", fit.delta)):=1:.N, by="id"]
      dt2[get(paste0("period", fit.delta))==1, (paste0("tstart", fit.delta)):=0]
      dt2[get(paste0("period", fit.delta))==1, (paste0("tstop", fit.delta)):=(get(time.var)<=fit.changepoint)*get(time.var)+
            (get(time.var)>fit.changepoint)*fit.changepoint]
      dt2[get(paste0("period", fit.delta))==1, (paste0("tstart", fit.delta)):=0]
      dt2[get(paste0("period", fit.delta))==1, (paste0("tstop", fit.delta)):=(get(time.var)<=fit.changepoint)*get(time.var)+
            (get(time.var)>fit.changepoint)*fit.changepoint]
      dt2[get(paste0("period", fit.delta))==2, (paste0("tstart", fit.delta)):=fit.changepoint]
      dt2[get(paste0("period", fit.delta))==2, (paste0("tstop", fit.delta)):=get(time.var)]
      dt2[get(paste0("period", fit.delta))==1 & !time.indicator, (delta.var):=delta1]
      mod1 <- as.character(fit.model)
      mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                 which(strsplit(mod1[2], "")[[1]]==",")-1), paste0("tstart", fit.delta, ", tstop", fit.delta), mod1[2]),
                     "~",
                     gsub(paste0("\\+", A.name, "\\+"), "", gsub(paste0("\\+", A.name, " "), "", gsub(" ", "", paste0("I((period", fit.delta,"==1)&(", A.name, "==1))",
                                                                                                                      " + I((period", fit.delta, "==2)&(", A.name, "==1))", " + ",
                                                                                                                      paste0("+",mod1[3]))))))
      fit.cox <- coxph(formula(mod2), data=dt2[!time.indicator | get(paste0("period", fit.delta))==1])
    } else {
      if (fit[1]=="sl" & length(grep("coxnet", sl.pick))>0) {
        X <- model.matrix(as.formula(deparse(fit.model)), data=dt)
        y <- dt[, Surv(get(time.var), get(delta.var)==fit.delta)]
        fit.cox <- glmnet(x=X, y=y, family="cox", maxit=1000,
                          lambda=fit.penalty)
      } else {
        fit.cox <- coxph(as.formula(paste0(deparse(fit.model), collapse="")), data=dt)
      }
    }

    estimation[[each]]$fit.cox <- fit.cox

    if (verbose) print(fit.cox)

    if (fit[1]=="km") {
      tmp <- suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE)))
      setnames(tmp, "strata", "A")
      tmp[, A:=as.numeric(gsub("A=", "", A))]
      bhaz.cox <- merge(bhaz.cox,
                        rbind(do.call("rbind", lapply(a, function(aa) data.table(time=0, hazard=0, A=aa))),
                              tmp),
                        by=c("time", "A"), all.x=TRUE)
      bhaz.cox[, hazard:=na.locf(hazard), by="A"]
      bhaz.cox[, (paste0("dhaz.", fit.delta)):=c(0, diff(hazard)), by="A"]
      setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
    } else {
      if (fit[1]=="sl" & length(grep("coxnet", sl.pick))>0) {
        basehaz <- glmnet_basesurv(dt[, get(time.var)],
                                   dt[, get(delta.var)==fit.delta], X, centered=TRUE)
        bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                          data.table(time=basehaz$time,
                                                     hazard=basehaz$cumulative_base_hazard)),
                          by="time", all.x=TRUE)
      } else {
        bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                          suppressWarnings(setDT(basehaz(fit.cox, centered=TRUE)))),
                          by="time", all.x=TRUE)
      }
      bhaz.cox[, hazard:=na.locf(hazard), by=A.name]
      bhaz.cox[, (paste0("dhaz", fit.delta)):=c(0, diff(hazard)), by=A.name]
      setnames(bhaz.cox, "hazard", paste0("chaz", fit.delta))
    }

  }

  setnames(bhaz.cox, c("time", A.name), c(time.var, A.name))

  bhaz.cox[, chaz0.1:=c(0, chaz0[-.N])]

  if (length(dt2)==0) {
    dt2 <- copy(dt)
  }

  any.hal <- unlist(lapply(estimation, function(each) length(grep("hal", each[["fit"]]))>0))

  if (!any(any.hal)) {
    bhaz.cox <- bhaz.cox[get(time.var)<=max(tau)]
  }

  dt2.a <- do.call("rbind", lapply(a, function(aa) {
    dt.tmp <- copy(dt2)
    dt.tmp[, (A.name):=aa]
  }))

  for (each in 1:length(estimation)) {

    if (estimation[[each]][["fit"]][1]=="sl" & length(grep("coxnet", estimation[[each]][["sl.pick"]]))>0) {
      X2.a <- model.matrix(as.formula(deparse(estimation[[each]][["model"]])), data=dt2.a)
      dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])):=
              exp(predict(estimation[[each]][["fit.cox"]], newx=X2.a, type="link"))]
    } else {
      dt2.a[, (paste0("fit.cox", estimation[[each]][["event"]])):=
              predict(estimation[[each]][["fit.cox"]], newdata=dt2.a, type="risk")]
    }

    dt2.a[get(paste0("fit.cox", estimation[[each]][["event"]]))>500, (paste0("fit.cox", estimation[[each]][["event"]])):=500]

  }

  mat <- do.call("rbind", lapply(1:n, function(i) {
    tmp <- cbind(id=i, bhaz.cox)
    tmp[, time.obs:=dt[id==i, get(time.var)]]
    tmp[, delta.obs:=dt[id==i, get(delta.var)]]
    tmp[, A.obs:=dt[id==i, get(A.name)]]
    tmp[, prob.A:=prob.A[i]]

    if (length(a)==2) {
      if (dt[, length(unique(get(A.name)))]>2) {
        tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) /
              ((prob.A))]
      } else {
        tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) /
              ((prob.A^get(A.name) * (1-prob.A)^(1-get(A.name))))]
      }
    } else {
      if (dt[, length(unique(get(A.name)))]>2) {
        tmp[, Ht:=-((get(A.name)==a)) /
              ((prob.A))]
      } else {
        tmp[, Ht:=-((get(A.name)==a)) /
              ((prob.A^get(A.name) * (1-prob.A)^(1-get(A.name))))]
      }
    }

    for (each in 1:length(estimation)) {
      if (length(estimation[[each]][["changepoint"]])>0) {
        tmp[, (paste0("period", estimation[[each]][["event"]])):=
              (get(time.var)<=estimation[[each]][["changepoint"]])*1+(get(time.var)>estimation[[each]][["changepoint"]])*2]
        tmp <- merge(tmp, dt2.a[id==i, c(paste0("period", estimation[[each]][["event"]]),
                                         A.name,
                                         paste0("fit.cox", estimation[[each]][["event"]])),
                                with=FALSE], by=c(paste0("period", estimation[[each]][["event"]]), A.name))
      } else {
        tmp <- merge(tmp, unique(dt2.a[id==i, c(A.name,
                                                paste0("fit.cox", estimation[[each]][["event"]])),
                                       with=FALSE]), by=c(A.name))
      }
    }

    if (any(unlist(lapply(estimation, function(x) x[["event"]]))==0)) {
      tmp[, surv.C1:=exp(-fit.cox0*chaz0.1)]
      if (tmp[, is.na(mean(is.na(surv.C1)))]) {
        stop("Censoring survival calculation returned all missing values.")
      }
      if (tmp[, any(is.na(surv.C1))]) {
        stop("Censoring survival calculation returned missing values.")
      }
      if (tmp[, min(surv.C1)<cut.off.cens]) {
        tmp[surv.C1<cut.off.cens, surv.C1:=cut.off.cens]
      }
      tmp[, Ht:=Ht/surv.C1]
    }
  }))

  if (verbose) paste0("min of censoring weights: ", mat[, min(surv.C1)])

  if (any(any.hal)) {

    set.seed(ceiling(sum(as.numeric(dt[[A.name]]))))

    count.hals <- 1

    mat <- merge(mat, dt[, c("id", c(covars[!covars%in%names(mat)])), with=FALSE], by="id")

    mat[, (covars):=lapply(.SD, function(x) {
      if (is.character(x)) {
        return(as.numeric(as.factor(x)))
      } else if (is.factor(x)) {
        return(as.numeric(x))
      } else return(x)
    }), .SDcols=covars]

    dt[, (covars):=lapply(.SD, function(x) {
      if (is.character(x)) {
        return(as.numeric(as.factor(x)))
      } else if (is.factor(x)) {
        return(as.numeric(x))
      } else return(x)
    }), .SDcols=covars]

    for (each in (1:length(estimation))[any.hal]) {

      covars1 <- covars

      fit.delta <- estimation[[each]][["event"]]
      fit.name <- names(estimation)[each]

      if (any(names(estimation[[each]])=="lambda.cvs")) {
        lambda.cvs.1 <- estimation[[each]][["lambda.cvs"]]
      } else {
        lambda.cvs.1 <- lambda.cvs
      }

      fit.type <- estimation[[each]][["fit"]][1]

      if (fit.type %in% c("cox.hal", "cox.hal.sl")) {

        if (do.hal.screening) {

          (one.way.screening <- hal.screening(covars=covars, dt=dt, cut.one.way=cut.one.way.screen,
                                              mat=mat, delta.var=delta.var, delta.value=fit.delta,
                                              treatment=A.name))

          (two.way.screening <- hal.screening(covars=one.way.screening, dt=dt, cut.one.way=cut.one.way.screen,
                                              delta.var=delta.var, delta.value=fit.delta,
                                              treatment=A.name, order=2,
                                              mat=mat, cut.two.way=cut.two.way.screen))

          if (verbose.hal) {
            print("one.way.screening:")
            print(one.way.screening)
            print("two.way.screening:")
            print(two.way.screening)
          }

          if (verbose & length(one.way.screening)>0) print(paste0("variables picked by initial screening: ", paste0(one.way.screening, collapse=", ")))
          mat <- fit.hal(covars=one.way.screening, dt=dt, cut.one.way=cut.one.way.final,
                         mat=mat, delta.var=delta.var, V=V,
                         intervention=a,
                         two.way=two.way.screening, cut.two.way=cut.two.way,
                         penalize.treatment=FALSE,
                         delta.value=fit.delta,
                         lambda.cvs=lambda.cvs.1,
                         verbose=verbose,
                         predict=max(tau), treatment.prediction=A.name,
                         treatment=A.name)

        } else  {
          mat <- fit.hal(covars=covars, dt=dt, cut.one.way=cut.one.way.final,
                         mat=mat, delta.var=delta.var, V=V,
                         intervention=a,
                         penalize.treatment=FALSE,
                         lambda.cvs=lambda.cvs.1,
                         verbose=verbose.hal, delta.value=fit.delta,
                         predict=max(tau), treatment.prediction=A.name,
                         treatment=A.name)

        }

      } else if (fit.type == "hal") {

        if (verbose) print(paste0("use HAL for = ", fit.name))

        if (do.hal.screening) {

          time.screen <- min(3, cut.time)
          one.way.screening <- hal.screening(covars=covars, dt=dt,
                                             cut.one.way=cut.one.way.screen,
                                             time.var=time.var, cut.time=time.screen,
                                             mat=mat, delta.var="delta.obs",
                                             delta.value=fit.delta,
                                             treatment="A.obs",
                                             cut.time.treatment=cut.time.A)

          two.way.screening <- cbind(var1="", var2="")
          if (length(one.way.screening)>0) {
            two.way.screening <- hal.screening(covars=one.way.screening, dt=dt,
                                               cut.one.way=cut.one.way.screen,
                                               time.var=time.var, cut.time=time.screen,
                                               delta.var="delta.obs",
                                               delta.value=fit.delta,
                                               treatment="A.obs",
                                               cut.time.treatment=cut.time.A, order=2,
                                               mat=mat, cut.two.way=cut.two.way.screen)
          }

          if (verbose.hal) {
            print("one.way.screening:")
            print(one.way.screening)
            print("two.way.screening:")
            print(two.way.screening)
          }
          if (verbose & length(one.way.screening)>0) print(paste0("variables picked by initial screening: ", paste0(one.way.screening, collapse=", ")))

          mat <- fit.hal(covars=one.way.screening, dt=dt, cut.one.way=cut.one.way.final,
                         time.var=time.var, cut.time=cut.time,
                         mat=mat, delta.var="delta.obs", V=V,
                         intervention=a,
                         verbose=verbose, delta.value=fit.delta,
                         two.way=two.way.screening, cut.two.way=cut.two.way,
                         penalize.treatment=FALSE, penalize.time=FALSE,
                         lambda.cvs=lambda.cvs.1,
                         predict=max(tau), treatment.prediction=A.name,
                         treatment="A.obs", cut.time.treatment=cut.time.A)

        } else if (hal.sl) {

          tune.grid <- expand.grid(cut.time=c(10, 13, 16),
                                   cut.time.treatment=c(0, 5, 10),
                                   cut.one.way=c(5, 10, 15),
                                   cut.two.way=c(0, 5, 10))

          cve.tune.grid <- sapply(1:nrow(tune.grid), function(jj) {
            return(fit.hal(covars=covars, dt=dt, cut.one.way=tune.grid[jj,"cut.one.way"],
                           time.var=time.var, cut.time=tune.grid[jj,"cut.time"],
                           mat=mat, delta.var="delta.obs", V=V,
                           intervention=a,
                           verbose=FALSE, delta.value=fit.delta,
                           cut.two.way=tune.grid[jj,"cut.two.way"],
                           penalize.treatment=FALSE, penalize.time=FALSE,
                           lambda.cvs=lambda.cvs.1,
                           predict=max(tau), treatment.prediction=A.name,
                           treatment="A.obs", cut.time.treatment=tune.grid[jj,"cut.time.treatment"],
                           return.cve=TRUE)$cve$min$cve)
          })

          tune.grid$cve <- cve.tune.grid
          tune.grid.min <- tune.grid[order(tune.grid$cut.time, tune.grid$cut.time.treatment,
                                           tune.grid$cut.one.way, tune.grid$cut.two.way), ]
          tune.grid.min <- tune.grid.min[tune.grid.min$cve==min(tune.grid.min$cve),][1,]

          mat <- fit.hal(covars=covars, dt=dt, cut.one.way=tune.grid.min[,"cut.one.way"],
                         time.var=time.var, cut.time=tune.grid.min[,"cut.time"],
                         mat=mat, delta.var="delta.obs", V=V,
                         intervention=a,
                         verbose=verbose.hal, delta.value=fit.delta,
                         cut.two.way=tune.grid.min[,"cut.two.way"],
                         penalize.treatment=FALSE, penalize.time=FALSE,
                         lambda.cvs=lambda.cvs.1,
                         predict=max(tau), treatment.prediction=A.name,
                         treatment="A.obs", cut.time.treatment=tune.grid.min[,"cut.time.treatment"])

        } else {

          mat <- fit.hal(covars=covars, dt=dt, cut.one.way=cut.one.way.final,
                         time.var=time.var, cut.time=cut.time,
                         mat=mat, delta.var="delta.obs", V=V,
                         intervention=a,
                         verbose=verbose.hal, delta.value=fit.delta,
                         penalize.treatment=FALSE, penalize.time=FALSE,
                         lambda.cvs=lambda.cvs.1,
                         predict=max(tau), treatment.prediction=A.name,
                         treatment="A.obs", cut.time.treatment=cut.time.A)

        }
      }

      count.hals <- count.hals+1
    }

    mat <- mat[get(time.var)<=max(tau)]

  }

  if (mat[, min(surv.C1)<cut.off.cens]) {
    if (mat[get(time.var)<=time.obs, min(surv.C1)<cut.off.cens]) {
      not.fit.list[[length(not.fit.list)+1]] <-
        paste0(paste0("there seems to be positivity issues for censoring; truncated ",
                      sum(mat[get(time.var)<=time.obs, sum(unique(surv.C1<cut.off.cens)), by="id"][,2][[1]]),
                      " observation(s), at level ",
                      cut.off.cens))
    }
    mat[surv.C1<=cut.off.cens, Ht:=Ht*surv.C1/cut.off.cens]
    mat[surv.C1<=cut.off.cens, surv.C1:=cut.off.cens]
  }

  if (mat[, min(prob.A)<cut.off.cens | max(prob.A)>(1-cut.off.cens)]) {
    if (mat[, min(prob.A)<cut.off.cens | max(prob.A)>(1-cut.off.cens)]) {
      not.fit.list[[length(not.fit.list)+1]] <-
        paste0(paste0("there seems to be positivity issues for treatment; truncated ",
                      sum(mat[, sum(unique(prob.A<cut.off.cens | prob.A>(1-cut.off.cens))), by="id"][,2][[1]]),
                      " observation(s), at level ",
                      cut.off.cens))
    }
    if (dt[, length(unique(get(A.name)))]>2) {
      mat[prob.A<cut.off.cens, Ht:=Ht*prob.A/cut.off.cens]
      mat[prob.A>(1-cut.off.cens), Ht:=Ht*prob.A/(1-cut.off.cens)]
    } else {
      mat[, prob.A.orig:=prob.A]
      mat[, prob.A.trunc:=pmin(pmax(prob.A.orig, cut.off.cens), 1-cut.off.cens)]
      mat[, denom.orig:=prob.A.orig^get(A.name) * (1-prob.A.orig)^(1-get(A.name))]
      mat[, denom.trunc:=prob.A.trunc^get(A.name) * (1-prob.A.trunc)^(1-get(A.name))]
      mat[denom.trunc>0, Ht:=Ht*denom.orig/denom.trunc]
      mat[, prob.A:=prob.A.trunc]
      mat[, c("prob.A.orig", "prob.A.trunc", "denom.orig", "denom.trunc"):=NULL]
    }
    if (dt[, length(unique(get(A.name)))]>2) {
      mat[prob.A<cut.off.cens, prob.A:=cut.off.cens]
      mat[prob.A>(1-cut.off.cens), prob.A:=1-cut.off.cens]
    }
  }

  mat[, surv.t:=1]
  for (each in outcome.index) {
    fit.delta <- estimation[[each]][["event"]]
    mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                       get(paste0("fit.cox", fit.delta)))),
        by=c("id", A.name)]
  }

  mat[, surv.t1:=c(1, surv.t[-.N]), by=c("id", A.name)]

  for (kk in 1:length(tau)) {
    mat[, (paste0("surv.tau", kk)):=
          surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
        by=c("id", A.name)]
  }

  for (each in outcome.index) {
    fit.delta <- estimation[[each]][["event"]]
    mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t1*get(paste0("dhaz", fit.delta))*
                                                   get(paste0("fit.cox", fit.delta))),
        by=c("id", A.name)]
    mat[, time.diff := get(time.var) - shift(get(time.var), fill=0), by=c("id", A.name)]
    mat[, (paste0("RMTL", fit.delta, ".t")) := cumsum(shift(get(paste0("F", fit.delta, ".t")), fill=0) * time.diff),
        by = c("id", A.name)]

    for (kk in 1:length(tau)) {
      mat[, (paste0("F", fit.delta, ".tau", kk)):=
            get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
          by=c("id", A.name)]
      mat[, (paste0("RMTL", fit.delta, ".tau", kk)) :=
            get(paste0("RMTL", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])]+
            get(paste0("F", fit.delta, ".tau",kk))*(tau[kk]-max(get(time.var)[get(time.var)<=tau[kk]])),
          by=c("id", A.name)]

      mat[, time_w := pmax(0, tau[kk] - get(time.var))]

      for (each2 in outcome.index) {
        fit.delta2 <- estimation[[each2]][["event"]]
        if (fit.delta==fit.delta2) {
          mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                -(time_w - (get(paste0("RMTL", fit.delta, ".tau", kk)) -
                              get(paste0("RMTL", fit.delta, ".t")) - get(paste0("F", fit.delta, ".t"))* time_w)  / surv.t) ]
          mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):= -time_w]
        } else {
          mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                (get(paste0("RMTL", fit.delta, ".tau", kk)) -
                   get(paste0("RMTL", fit.delta, ".t")) - get(paste0("F", fit.delta, ".t"))* time_w) / surv.t ]
          mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=0]
        }
      }
    }

  }

  init.fit.rmtl <- lapply(target, function(each) {
    sapply(1:length(tau), function(kk) {

      mean(rowSums(sapply(a, function(aa)
        (2*(aa==a[1])-1) * (
          mat[get(A.name)==aa,
              get(paste0("RMTL", estimation[[outcome.index[each]]][["event"]],
                         ".tau", kk))[1],
              by="id"][,2][[1]]
        )
      )))
    })
  })
  names(init.fit.rmtl) <- paste0("RMTL", sapply(outcome.index[target],
                                                function(each) estimation[[each]][["event"]]))

  eval.ic.rmtl <- function(mat, fit, target.index=outcome.index, Sigma=FALSE, tau.values=tau,
                           survival=FALSE, tau.all=tau, return.ic=FALSE) {
    time.var <- trimws(as.character(time.var)[1L])
    tau.values <- as.numeric(unlist(tau.values))
    tau.all <- as.numeric(unlist(tau.all))
    if (!time.var %in% names(mat)) stop("Time variable not found in mat: ", time.var)
    if (!is.atomic(mat[[time.var]])) {
      stop("Time variable is not atomic: ", time.var,
           "; class = ", paste(class(mat[[time.var]]), collapse=", "))
    }
    outer <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      each.index <- (1:length(target.index))[target.index==each]
      ic.mat <- sapply(1:length(tau.values), function(kk) {
        k2 <- (1:length(tau.all))[tau.all==max(tau.all[tau.all<=tau.values[kk]])]
        target_tau <- tau.values[kk]

        out_rmtl <- 0
        for (each2 in outcome.index) {

          fit.delta2 <- estimation[[each2]][["event"]]
          mat[get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2))<=-500,
              (paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)):=-500]

          out_rmtl2 <- mat[, {
            time_col <- .SD[[time.var]]
            sum( (.SD[[A.name]]==A.obs) * (time_col<=tau.values[kk]) *
                   (time_col<=time.obs) * Ht *
                   get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)) *
                   ( (delta.obs==fit.delta2 & time_col==time.obs) -
                       get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2)) ),
                 na.rm=TRUE)
          }, by="id"]
          out_rmtl <- out_rmtl +  out_rmtl2[, 2][[1]]
        }

        ic.squared <- (out_rmtl + rowSums(sapply(a, function(aa)
          (2*(aa==a[1])-1)*(mat[get(A.name)==aa, get(paste0("RMTL", fit.delta,
                                                            ".tau", k2))[1],
                                by="id"][,2][[1]])))-
            (length(a)==2)*fit[[paste0("RMTL", fit.delta)]][kk]-
            (length(a)==1)*fit[[paste0("RMTL", fit.delta)]][kk])
        if (survival | Sigma | return.ic) return(ic.squared) else return(sqrt(mean(ic.squared^2)/n))
      })
      if ((survival | Sigma | return.ic) & is.null(dim(ic.mat))) {
        ic.mat <- matrix(ic.mat, ncol=1)
      }
      if (survival | Sigma | return.ic) {
        colnames(ic.mat) <- paste0("tau=", tau.values)
      }
      return(ic.mat)
    })
    if (return.ic) {
      names(outer) <- paste0("RMTL", sapply(target.index, function(each) estimation[[each]]["event"]))
      return(outer)
    }
    if (Sigma) {
      if (survival & length(outcome.index)==length(target)) {
        outer2 <- Reduce("+", lapply(outer, function(xx) -xx))
      } else {
        outer2 <- do.call("cbind", outer)
      }
      Sigma.list <- lapply(1:n, function(i) {
        t(outer2[i,,drop=FALSE])%*%outer2[i,,drop=FALSE]
      })
      return(Reduce("+", Sigma.list) / length(Sigma.list))
    } else {
      names(outer) <- paste0("RMTL", sapply(target.index, function(each) estimation[[each]]["event"]))
      return(outer)
    }
  }
  combine.rmst.free.ic <- function(ic.list) {
    mats <- lapply(ic.list, function(x) {
      if (is.null(dim(x))) matrix(x, ncol=1) else as.matrix(x)
    })
    out <- Reduce("+", lapply(mats, function(x) -x))
    colnames(out) <- colnames(mats[[1]])
    out
  }
  init.ic.rmtl <- eval.ic.rmtl(mat, fit=init.fit.rmtl, target.index=outcome.index[target])
  init.ic.rmtl.full <- eval.ic.rmtl(mat, fit=init.fit.rmtl,
                                    target.index=outcome.index[target],
                                    return.ic=TRUE)
  init.se.rmst.free <- NULL
  if (length(outcome.index)>1 && length(target)==length(outcome.index) && !target.S) {
    init.ic.rmst.free <- combine.rmst.free.ic(init.ic.rmtl.full)
    init.se.rmst.free <- sqrt(colMeans(init.ic.rmst.free^2)/n)
  }

  if (weighted.norm=="Sigma") {
    Sigma <- eval.ic.rmtl(mat, fit=init.fit.rmtl, target.index=outcome.index[target], Sigma=TRUE)
    Sigma.inv <- try(solve(Sigma), silent=TRUE)
    if (inherits(Sigma.inv, "try-error")) {
      Sigma.inv <- solve(Sigma + diag(x=0.000001, nrow=nrow(Sigma)))
      warning("regularization of Sigma needed for inversion")
      not.fit.list[[length(not.fit.list)+1]] <- paste0("Sigma was regularized")
    }
  }

  init.list <-  lapply(1:length(init.fit.rmtl), function(each.index) {
    out <- rbind(init.est=init.fit.rmtl[[each.index]],
                 init.se=init.ic.rmtl[[each.index]])
    colnames(out) <- paste0("tau=", tau)
    return(out)
  })
  names(init.list) <- paste0("RMTL", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))

  tmle.list <- list(init=init.list)

  summarize.arm.rmtl.fit <- function(mat, target.index=outcome.index[target],
                                     tau.values=tau, tau.all=tau, arm.values=a) {
    out <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      fit.mat <- do.call("rbind", lapply(arm.values, function(aa) {
        sapply(seq_along(tau.values), function(kk) {
          k2 <- max(which(tau.all <= tau.values[kk]))
          mean(mat[get(A.name)==aa,
                   get(paste0("RMTL", fit.delta, ".tau", k2))[1],
                   by="id"][, 2][[1]])
        })
      }))
      rownames(fit.mat) <- as.character(arm.values)
      colnames(fit.mat) <- paste0("tau=", tau.values)
      fit.mat
    })
    names(out) <- paste0("RMTL", sapply(target.index, function(each) estimation[[each]][["event"]]))
    out
  }

  summarize.arm.cif.fit <- function(mat, target.index=outcome.index[target],
                                    tau.values=tau, tau.all=tau, arm.values=a) {
    out <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      fit.mat <- do.call("rbind", lapply(arm.values, function(aa) {
        sapply(seq_along(tau.values), function(kk) {
          k2 <- max(which(tau.all <= tau.values[kk]))
          mean(mat[get(A.name)==aa,
                   get(paste0("F", fit.delta, ".tau", k2))[1],
                   by="id"][, 2][[1]])
        })
      }))
      rownames(fit.mat) <- as.character(arm.values)
      colnames(fit.mat) <- paste0("tau=", tau.values)
      fit.mat
    })
    names(out) <- paste0("F", sapply(target.index, function(each) estimation[[each]][["event"]]))
    out
  }

  eval.ic.rmtl.arm <- function(mat, fit, target.index=outcome.index[target],
                               tau.values=tau, tau.all=tau, arm.values=a,
                               return.ic=FALSE) {
    time.var <- trimws(as.character(time.var)[1L])
    tau.values <- as.numeric(unlist(tau.values))
    tau.all <- as.numeric(unlist(tau.all))
    if (!time.var %in% names(mat)) stop("Time variable not found in mat: ", time.var)
    if (!is.atomic(mat[[time.var]])) {
      stop("Time variable is not atomic: ", time.var,
           "; class = ", paste(class(mat[[time.var]]), collapse=", "))
    }
    out <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      ic.by.arm <- lapply(seq_along(arm.values), function(aa.index) {
        aa <- arm.values[aa.index]
        ic.mat <- sapply(seq_along(tau.values), function(kk) {
          k2 <- max(which(tau.all <= tau.values[kk]))
          plug.in <- mat[get(A.name)==aa,
                         get(paste0("RMTL", fit.delta, ".tau", k2))[1],
                         by="id"][order(id), 2][[1]]
          ic.correction <- numeric(n)

          for (each2 in outcome.index) {
            fit.delta2 <- estimation[[each2]][["event"]]
            tmp <- mat[mat[[A.name]]==aa & mat[["A.obs"]]==aa &
                         mat[[time.var]]<=tau.values[kk] &
                         mat[[time.var]]<=mat[["time.obs"]],
                       {
                         time_col <- .SD[[time.var]]
                         .(term=sum(abs(Ht) *
                                      get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", k2)) *
                                      (((delta.obs==fit.delta2) & (time_col==time.obs)) -
                                         get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2))),
                                    na.rm=TRUE))
                       },
                       by="id"]
            if (nrow(tmp)>0) {
              ic.correction[tmp$id] <- ic.correction[tmp$id] + tmp$term
            }
          }

          ic.correction + plug.in - fit[[paste0("RMTL", fit.delta)]][aa.index, kk]
        })
        if (is.null(dim(ic.mat))) {
          ic.mat <- matrix(ic.mat, ncol=1)
        }
        colnames(ic.mat) <- paste0("tau=", tau.values)
        ic.mat
      })
      names(ic.by.arm) <- as.character(arm.values)

      if (return.ic) {
        return(ic.by.arm)
      }

      se.mat <- do.call("rbind", lapply(ic.by.arm, function(ic.mat) {
        sqrt(colMeans(ic.mat^2)/n)
      }))
      rownames(se.mat) <- as.character(arm.values)
      colnames(se.mat) <- paste0("tau=", tau.values)
      se.mat
    })
    names(out) <- paste0("RMTL", sapply(target.index, function(each) estimation[[each]][["event"]]))
    out
  }

  summarize.arm.rmst.free.fit <- function(fit, target.index=outcome.index[target],
                                          tau.values=tau, arm.values=a) {
    fit.mat <- do.call("rbind", lapply(seq_along(arm.values), function(aa.index) {
      tau.values - Reduce("+", lapply(target.index, function(each) {
        fit[[paste0("RMTL", estimation[[each]][["event"]])]][aa.index, ]
      }))
    }))
    rownames(fit.mat) <- as.character(arm.values)
    colnames(fit.mat) <- paste0("tau=", tau.values)
    list("RMST.free"=fit.mat)
  }

  combine.arm.rmst.free.ic <- function(ic.by.event, arm.values=a) {
    out <- lapply(as.character(arm.values), function(aa) {
      combine.rmst.free.ic(lapply(ic.by.event, function(x) x[[aa]]))
    })
    names(out) <- as.character(arm.values)
    out
  }

  summarize.arm.se.from.ic <- function(ic.by.event, arm.values=a) {
    out <- lapply(names(ic.by.event), function(name.now) {
      se.mat <- do.call("rbind", lapply(ic.by.event[[name.now]], function(ic.mat) {
        sqrt(colMeans(ic.mat^2)/n)
      }))
      rownames(se.mat) <- as.character(arm.values)
      colnames(se.mat) <- colnames(ic.by.event[[name.now]][[1]])
      se.mat
    })
    names(out) <- names(ic.by.event)
    out
  }

  eval.ic.cif.arm <- function(mat, fit, target.index=outcome.index[target],
                              tau.values=tau, tau.all=tau, arm.values=a,
                              return.ic=FALSE) {
    time.var <- trimws(as.character(time.var)[1L])
    tau.values <- as.numeric(unlist(tau.values))
    tau.all <- as.numeric(unlist(tau.all))
    if (!time.var %in% names(mat)) stop("Time variable not found in mat: ", time.var)
    if (!is.atomic(mat[[time.var]])) {
      stop("Time variable is not atomic: ", time.var,
           "; class = ", paste(class(mat[[time.var]]), collapse=", "))
    }
    out <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      ic.by.arm <- lapply(seq_along(arm.values), function(aa.index) {
        aa <- arm.values[aa.index]
        ic.mat <- sapply(seq_along(tau.values), function(kk) {
          k2 <- max(which(tau.all <= tau.values[kk]))
          plug.in <- mat[get(A.name)==aa,
                         get(paste0("F", fit.delta, ".tau", k2))[1],
                         by="id"][order(id), 2][[1]]
          ic.correction <- numeric(n)

          for (each2 in outcome.index) {
            fit.delta2 <- estimation[[each2]][["event"]]
            tmp <- mat[mat[[A.name]]==aa & mat[["A.obs"]]==aa &
                         mat[[time.var]]<=tau.values[kk] &
                         mat[[time.var]]<=mat[["time.obs"]],
                       {
                         time_col <- .SD[[time.var]]
                         .(term=sum(abs(Ht) *
                                      get(paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", k2)) *
                                      (((delta.obs==fit.delta2) & (time_col==time.obs)) -
                                         get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2))),
                                    na.rm=TRUE))
                       },
                       by="id"]
            if (nrow(tmp)>0) {
              ic.correction[tmp$id] <- ic.correction[tmp$id] + tmp$term
            }
          }

          ic.correction + plug.in - fit[[paste0("F", fit.delta)]][aa.index, kk]
        })
        if (is.null(dim(ic.mat))) {
          ic.mat <- matrix(ic.mat, ncol=1)
        }
        colnames(ic.mat) <- paste0("tau=", tau.values)
        ic.mat
      })
      names(ic.by.arm) <- as.character(arm.values)

      if (return.ic) {
        return(ic.by.arm)
      }

      se.mat <- do.call("rbind", lapply(ic.by.arm, function(ic.mat) {
        sqrt(colMeans(ic.mat^2)/n)
      }))
      rownames(se.mat) <- as.character(arm.values)
      colnames(se.mat) <- paste0("tau=", tau.values)
      se.mat
    })
    names(out) <- paste0("F", sapply(target.index, function(each) estimation[[each]][["event"]]))
    out
  }

  simulate.qmax <- function(rho, repeat.no=50000, seed=1) {
    if (nrow(rho)==1) {
      return(qnorm(0.975))
    }
    rho[is.na(rho)] <- 0
    diag(rho) <- 1
    set.seed(seed)
    z.sim <- try(MASS::mvrnorm(n=repeat.no, mu=rep(0, nrow(rho)), Sigma=rho,
                               tol=1e-6, empirical=FALSE), silent=TRUE)
    if (inherits(z.sim, "try-error")) {
      warning("Could not simulate simultaneous critical value; returning NA")
      return(NA_real_)
    }
    as.numeric(stats::quantile(apply(abs(z.sim), 1, max), probs=0.95, na.rm=TRUE))
  }

  compute.simultaneous.cif.band <- function(ic, fit, target.index=outcome.index[target],
                                            tau.values=tau, arm.values=a) {
    out <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      q.vec <- setNames(rep(NA_real_, length(arm.values)), as.character(arm.values))
      lower.mat <- upper.mat <- matrix(NA_real_, nrow=length(arm.values), ncol=length(tau.values))
      rownames(lower.mat) <- rownames(upper.mat) <- as.character(arm.values)
      colnames(lower.mat) <- colnames(upper.mat) <- paste0("tau=", tau.values)

      for (aa.index in seq_along(arm.values)) {
        aa <- arm.values[aa.index]
        ic.mat <- ic[[paste0("F", fit.delta)]][[as.character(aa)]]
        sigma <- crossprod(ic.mat) / nrow(ic.mat)
        sigma.diag <- diag(sigma)
        if (all(sigma.diag == 0)) {
          q.now <- qnorm(0.975)
        } else if (any(is.na(sigma.diag)) || any(sigma.diag < 0) || any(sigma.diag == 0)) {
          q.now <- NA_real_
        } else {
          denom <- sqrt(outer(sigma.diag, sigma.diag))
          rho <- sigma / denom
          q.now <- simulate.qmax(rho,
                                 seed=ceiling(sum(as.numeric(dt[[A.name]]))) +
                                   100*fit.delta + aa.index)
        }
        q.vec[as.character(aa)] <- q.now
        se.vec <- sqrt(diag(sigma)/n)
        lower.mat[aa.index, ] <- fit[[paste0("F", fit.delta)]][aa.index, ] - q.now*se.vec
        upper.mat[aa.index, ] <- fit[[paste0("F", fit.delta)]][aa.index, ] + q.now*se.vec
      }

      list(q.max.95=q.vec, lower=lower.mat, upper=upper.mat)
    })
    names(out) <- paste0("F", sapply(target.index, function(each) estimation[[each]][["event"]]))
    out
  }

  compute.simultaneous.rmtld.band <- function(ic, fit, target.index=outcome.index[target],
                                              tau.values=tau) {
    out <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      ic.mat <- ic[[paste0("RMTL", fit.delta)]]
      sigma <- crossprod(ic.mat) / nrow(ic.mat)
      sigma.diag <- diag(sigma)
      if (all(sigma.diag == 0)) {
        q.now <- qnorm(0.975)
      } else if (any(is.na(sigma.diag)) || any(sigma.diag < 0) || any(sigma.diag == 0)) {
        q.now <- NA_real_
      } else {
        denom <- sqrt(outer(sigma.diag, sigma.diag))
        rho <- sigma / denom
        q.now <- simulate.qmax(rho,
                               seed=ceiling(sum(as.numeric(dt[[A.name]]))) +
                                 1000 + fit.delta)
      }
      se.vec <- sqrt(diag(sigma)/n)
      list(q.max.95=q.now,
           lower=fit[[paste0("RMTL", fit.delta)]] - q.now*se.vec,
           upper=fit[[paste0("RMTL", fit.delta)]] + q.now*se.vec)
    })
    names(out) <- paste0("RMTL", sapply(target.index, function(each) estimation[[each]][["event"]]))
    out
  }

  compute.simultaneous.scalar.band <- function(ic.mat, fit.vec, seed.offset=0) {
    sigma <- crossprod(ic.mat) / nrow(ic.mat)
    sigma.diag <- diag(sigma)
    if (all(sigma.diag == 0)) {
      q.now <- qnorm(0.975)
    } else if (any(is.na(sigma.diag)) || any(sigma.diag < 0) || any(sigma.diag == 0)) {
      q.now <- NA_real_
    } else {
      denom <- sqrt(outer(sigma.diag, sigma.diag))
      rho <- sigma / denom
      q.now <- simulate.qmax(rho,
                             seed=ceiling(sum(as.numeric(dt[[A.name]]))) + seed.offset)
    }
    se.vec <- sqrt(diag(sigma)/n)
    list(q.max.95=q.now,
         lower=fit.vec - q.now*se.vec,
         upper=fit.vec + q.now*se.vec)
  }

  format.arm.rmtl.output <- function(fit, se, target.index=outcome.index[target],
                                     tau.values=tau, arm.values=a) {
    data.table::rbindlist(lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      out <- expand.grid(arm=arm.values, tau=tau.values, stringsAsFactors=FALSE)
      out$event <- fit.delta
      out$est <- as.vector(fit[[paste0("RMTL", fit.delta)]])
      out$se <- as.vector(se[[paste0("RMTL", fit.delta)]])
      out$lower <- out$est - 1.96*out$se
      out$upper <- out$est + 1.96*out$se
      data.table(out)[, .(event, arm, tau, est, se, lower, upper)]
    }))
  }

  format.tmle.rmtld.output <- function(fit, se, sim.band=NULL,
                                       target.index=outcome.index[target],
                                       tau.values=tau) {
    data.table::rbindlist(lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      out <- data.frame(event=fit.delta, tau=tau.values, stringsAsFactors=FALSE)
      out$est <- as.vector(fit[[paste0("RMTL", fit.delta)]])
      out$se <- as.vector(se[[paste0("RMTL", fit.delta)]])
      out$lower <- out$est - 1.96*out$se
      out$upper <- out$est + 1.96*out$se
      out$sim.lower <- NA_real_
      out$sim.upper <- NA_real_
      out$q.max.95 <- NA_real_

      if (!is.null(sim.band)) {
        sim.now <- sim.band[[paste0("RMTL", fit.delta)]]
        if (!is.null(sim.now)) {
          out$sim.lower <- as.vector(sim.now$lower)
          out$sim.upper <- as.vector(sim.now$upper)
          out$q.max.95 <- sim.now$q.max.95
        }
      }

      data.table(out)[, .(event, tau, est, se, lower, upper,
                          sim.lower, sim.upper, q.max.95)]
    }))
  }

  format.tmle.rmst.free.output <- function(fit, se, sim.band=NULL,
                                           tau.values=tau) {
    out <- data.frame(tau=tau.values, stringsAsFactors=FALSE)
    out$est <- as.vector(fit[["RMST.free"]])
    out$se <- as.vector(se[["RMST.free"]])
    out$lower <- out$est - 1.96*out$se
    out$upper <- out$est + 1.96*out$se
    out$sim.lower <- NA_real_
    out$sim.upper <- NA_real_
    out$q.max.95 <- NA_real_
    if (!is.null(sim.band)) {
      out$sim.lower <- as.vector(sim.band$lower)
      out$sim.upper <- as.vector(sim.band$upper)
      out$q.max.95 <- sim.band$q.max.95
    }
    data.table(out)[, .(tau, est, se, lower, upper, sim.lower, sim.upper, q.max.95)]
  }

  format.arm.rmst.free.output <- function(fit, se, tau.values=tau, arm.values=a) {
    out <- expand.grid(arm=arm.values, tau=tau.values, stringsAsFactors=FALSE)
    out$est <- as.vector(fit[["RMST.free"]])
    out$se <- as.vector(se[["RMST.free"]])
    out$lower <- out$est - 1.96*out$se
    out$upper <- out$est + 1.96*out$se
    data.table(out)[, .(arm, tau, est, se, lower, upper)]
  }

  format.arm.cif.output <- function(fit, se,
                                    target.index=outcome.index[target],
                                    tau.values=tau, arm.values=a) {
    data.table::rbindlist(lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      out <- expand.grid(arm=arm.values, tau=tau.values, stringsAsFactors=FALSE)
      out$event <- fit.delta
      out$est <- as.vector(fit[[paste0("F", fit.delta)]])
      out$se <- as.vector(se[[paste0("F", fit.delta)]])
      out$lower <- out$est - 1.96*out$se
      out$upper <- out$est + 1.96*out$se
      data.table(out)[, .(event, arm, tau, est, se, lower, upper)]
    }))
  }

  summarize.adjusted.cif <- function(mat, target.index=outcome.index[target], arm.values=a) {
    data.table::rbindlist(lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      tmp <- mat[get(A.name) %in% arm.values,
                 .(cif=mean(get(paste0("F", fit.delta, ".t"))),
                   rmtl=mean(get(paste0("RMTL", fit.delta, ".t")))),
                 by=.(arm=get(A.name), time=get(time.var))]
      tmp[, event:=fit.delta]
      tmp[, .(event, arm, time, cif, rmtl)][order(event, arm, time)]
    }))
  }

  eval.equation <- function(mat, eps=0, target.index=outcome.index, cr.index=outcome.index,
                            tau.values=tau, tau.all=tau) {
    time.var <- trimws(as.character(time.var)[1L])
    tau.values <- as.numeric(unlist(tau.values))
    tau.all <- as.numeric(unlist(tau.all))
    if (!time.var %in% names(mat)) stop("Time variable not found in mat: ", time.var)
    if (!is.atomic(mat[[time.var]])) {
      stop("Time variable is not atomic: ", time.var,
           "; class = ", paste(class(mat[[time.var]]), collapse=", "))
    }
    outer <- lapply(target.index, function(each) {
      fit.delta <- estimation[[each]][["event"]]
      sapply(1:length(tau.values), function(kk) {
        k2 <- (1:length(tau.all))[tau.all==max(tau.all[tau.all<=tau.values[kk]])]
        out <- 0
        for (each2 in cr.index) {
          fit.delta2 <- estimation[[each2]][["event"]]
          out2 <- mat[(mat[[time.var]]<=tau.values[kk]) & (mat[[time.var]]<=mat[["time.obs"]]),
                      {
                        time_col <- .SD[[time.var]]
                        sum( (.SD[[A.name]]==A.obs) *
                               (time_col<=time.obs) * Ht *
                               get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)) *
                               ( (delta.obs==fit.delta2 & time_col==time.obs) -
                                   exp( eps * Ht *
                                          get(paste0("Ht", fit.delta, ".lambda", fit.delta2,".", k2)) ) *
                                   get(paste0("dhaz", fit.delta2)) * get(paste0("fit.cox", fit.delta2)) ))
                      }, by="id"]
          out <- out + out2[, 2][[1]]
        }
        return(mean(out))
      })
    })
    names(outer) <- paste0("RMTL", sapply(target.index, function(each) estimation[[each]]["event"]))
    return(outer)
  }

  Pn.eic.fun <- function(mat) {
    eval <- eval.equation(mat, eps=0, target.index=outcome.index[target])
    if (target.S) {
      eval <- sapply(tau, function(tt) -sum(sapply(eval, function(xx) xx[tau==tt])))
    }
    return(eval)
  }

  current.fit.rmtl.fun <- function(mat) {
    current.fit <- lapply(outcome.index[target], function(each) {
      sapply(1:length(tau), function(kk) {
        mean(rowSums(sapply(a, function(aa)
          (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                get(paste0("RMTL", estimation[[each]][["event"]],
                                           ".tau", kk))[1],
                                by="id"][,2][[1]]))))
      })
    })
    names(current.fit) <- paste0("RMTL", sapply(outcome.index[target], function(each) estimation[[each]][["event"]]))
    current.fit
  }

  eic.convergence.scale.fun <- function(mat) {
    if (eic.scale=="initial") {
      return(list(rmtl=init.ic.rmtl, rmst.free=init.se.rmst.free))
    }
    current.fit <- current.fit.rmtl.fun(mat)
    current.ic.rmtl <- eval.ic.rmtl(mat, fit=current.fit, target.index=outcome.index[target])
    current.se.rmst.free <- NULL
    if (length(outcome.index)>1 && length(target)==length(outcome.index) && !target.S) {
      current.ic.rmtl.full <- eval.ic.rmtl(mat, fit=current.fit,
                                           target.index=outcome.index[target],
                                           return.ic=TRUE)
      current.se.rmst.free <- sqrt(colMeans(combine.rmst.free.ic(current.ic.rmtl.full)^2)/n)
    }
    list(rmtl=current.ic.rmtl, rmst.free=current.se.rmst.free)
  }

  standardize.eic.fun <- function(Pn.eic, scale) {
    lapply(1:length(Pn.eic), function(kk) {
      denom <- scale[[kk]]
      if (any(unlist(scale)==0)) denom <- denom + 0.001
      Pn.eic[[kk]] / (denom * sqrt(n))
    })
  }

  max.standardized.eic.fun <- function(Pn.eic, scale.list) {
    Pn.eic.std <- standardize.eic.fun(Pn.eic, scale.list$rmtl)
    if (length(outcome.index)>1 && length(target)==length(outcome.index) &&
        !target.S && !is.null(scale.list$rmst.free)) {
      Pn.eic.std[[length(Pn.eic.std)+1]] <-
        -Reduce("+", Pn.eic) /
        ((ifelse(any(scale.list$rmst.free==0),
                 scale.list$rmst.free+0.001,
                 scale.list$rmst.free)) * sqrt(n))
    }
    max(abs(unlist(Pn.eic.std)))
  }

  Pn.eic <- Pn.eic.fun(mat)

  Pn.eic.norm.fun <- function(x2, x) {
    return(sqrt(sum(unlist(x2)*unlist(x))))
  }

  if (weighted.norm=="sigma") {
    Pn.eic2.fun <- function(Pn.eic) {
      if (is.list(Pn.eic)) {
        lapply(1:length(Pn.eic), function(kk)
          Pn.eic[[kk]] / ((ifelse(any(unlist(init.ic.rmtl)==0), init.ic.rmtl[[kk]]+0.00001, init.ic.rmtl[[kk]])^2) * n))
      } else {
        Pn.eic / ((ifelse(any(unlist(init.ic.rmtl)==0), init.ic.rmtl+0.00001, init.ic.rmtl)^2) * n)
      }
    }
  } else if (weighted.norm=="Sigma") {
    Pn.eic2.fun <- function(Pn.eic) {
      Pn.eic.tmp <- as.numeric(matrix(unlist(Pn.eic), 1, length(unlist(Pn.eic))) %*% Sigma.inv)
      if (length(target)>1) {
        lapply(1:length(target), function(kk) Pn.eic.tmp[(1+(kk-1)*length(tau)):(kk*length(tau))])
      } else {
        Pn.eic.tmp
      }
    }
  } else {
    Pn.eic2.fun <- function(Pn.eic) Pn.eic
  }

  Pn.eic2 <- Pn.eic2.fun(Pn.eic)

  criterion <- 1/(sqrt(n)*log(n))

  Pn.eic.norm.prev <- Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)

  if (verbose) print(paste0("Pn.eic.norm = ", Pn.eic.norm))

  for (step in 1:no.small.steps) {

    for (each in outcome.index) {
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("fit.cox", fit.delta, ".tmp")):=get(paste0("fit.cox", fit.delta))]
      for (kk in 1:length(tau)) {
        if (target.S) mat[, (paste0("Ht", ".lambda", fit.delta,".", kk, ".tmp")):=
                            get(paste0("Ht", ".lambda", fit.delta,".", kk))]
        for (each2 in outcome.index[target]) {
          fit.delta2 <- estimation[[each2]][["event"]]
          mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk, ".tmp")):=
                get(paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk))]
        }
      }
    }

    for (each in outcome.index) {
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("delta", fit.delta, ".dx")):=0]
      for (each2 in outcome.index[target]) {
        fit.delta2 <- estimation[[each2]][["event"]]
        for (kk in 1:length(tau)) {
          mat[, (paste0("delta", fit.delta, ".dx")):=
                get(paste0("delta", fit.delta, ".dx"))+
                (get(time.var)<=tau[kk])*Ht*(
                  (get(paste0("Ht", fit.delta2,".lambda", fit.delta,".", kk)))*
                    unlist(Pn.eic2[each2==outcome.index[target]])[kk]
                )/Pn.eic.norm]
        }
      }
    }

    deps <- deps.size

    mat[, surv.t:=1]
    for (each in outcome.index) {
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("fit.cox", fit.delta)):=
            get(paste0("fit.cox", fit.delta))*
            exp(deps*get(paste0("delta", fit.delta, ".dx")))]
      mat[get(paste0("fit.cox", fit.delta))>500, (paste0("fit.cox", fit.delta)):=500]
      mat[, surv.t:=surv.t*exp(-cumsum(get(paste0("dhaz", fit.delta))*
                                         get(paste0("fit.cox", fit.delta)))),
          by=c("id", A.name)]
      mat[get(paste0("fit.cox", fit.delta))==Inf, surv.t:=0]
      mat[surv.t < 1e-8,surv.t:=0]
    }

    mat[, surv.t1:=c(1, surv.t[-.N]), by=c("id", A.name)]
    for (kk in 1:length(tau)) {
      mat[, (paste0("surv.tau", kk)):=
            surv.t[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
          by=c("id", A.name)]
    }

    for (each in outcome.index[target]) {
      fit.delta <- estimation[[each]][["event"]]
      mat[, (paste0("F", fit.delta, ".t")):=cumsum(surv.t1*get(paste0("dhaz", fit.delta))*
                                                     get(paste0("fit.cox", fit.delta))),
          by=c("id", A.name)]
      mat[, time.diff := get(time.var) - shift(get(time.var), fill=0), by=c("id", A.name)]
      mat[, (paste0("RMTL", fit.delta, ".t")) := cumsum(shift(get(paste0("F", fit.delta, ".t")), fill=0) * time.diff),
          by = c("id", A.name)]

      for (kk in 1:length(tau)) {
        mat[, (paste0("F", fit.delta, ".tau", kk)):=
              get(paste0("F", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])],
            by=c("id", A.name)]
        mat[, (paste0("RMTL", fit.delta, ".tau", kk)) :=
              get(paste0("RMTL", fit.delta, ".t"))[get(time.var)==max(get(time.var)[get(time.var)<=tau[kk]])]+
              get(paste0("F", fit.delta, ".tau",kk))*(tau[kk]-max(get(time.var)[get(time.var)<=tau[kk]])),
            by=c("id", A.name)]

        mat[, time_w := pmax(0, tau[kk] - get(time.var))]

        for (each2 in outcome.index) {
          fit.delta2 <- estimation[[each2]][["event"]]
          if (fit.delta==fit.delta2) {
            mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                  -(time_w - (get(paste0("RMTL", fit.delta, ".tau", kk)) -
                                get(paste0("RMTL", fit.delta, ".t")) - get(paste0("F", fit.delta, ".t"))* time_w)  / surv.t) ]
            mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):= -time_w]
          } else {
            mat[surv.t>0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=
                  (get(paste0("RMTL", fit.delta, ".tau", kk)) -
                     get(paste0("RMTL", fit.delta, ".t")) - get(paste0("F", fit.delta, ".t"))* time_w) / surv.t ]
            mat[surv.t==0, (paste0("Ht", fit.delta, ".lambda", fit.delta2, ".", kk)):=0]
          }
        }
      }
    }

    Pn.eic <- Pn.eic.fun(mat)

    Pn.eic2 <- Pn.eic2.fun(Pn.eic)

    Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)

    if (verbose) print(paste0("step = ", step))

    if (verbose) print(paste0("Pn.eic.norm=", Pn.eic.norm))

    if (Pn.eic.norm.prev<=Pn.eic.norm) {

      if (verbose) {
        print("----")
        print(step)
        print("----")
      }

      for (each in outcome.index) {
        fit.delta <- estimation[[each]][["event"]]
        mat[, (paste0("fit.cox", fit.delta)):=get(paste0("fit.cox", fit.delta, ".tmp"))]
        for (kk in 1:length(tau)) {
          if (target.S) mat[, (paste0("Ht", ".lambda", fit.delta,".", kk)):=
                              get(paste0("Ht", ".lambda", fit.delta,".", kk, ".tmp"))]
          for (each2 in outcome.index[target]) {
            fit.delta2 <- estimation[[each2]][["event"]]
            mat[, (paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk)):=
                  get(paste0("Ht", fit.delta2, ".lambda", fit.delta,".", kk, ".tmp"))]
          }
        }
      }

      Pn.eic <- Pn.eic.fun(mat)
      Pn.eic2 <- Pn.eic2.fun(Pn.eic)

      Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
      deps.size <- 0.5*deps.size

    } else {

      Pn.eic.norm.prev <- Pn.eic.norm

    }

    eic.convergence.scale <- eic.convergence.scale.fun(mat)
    Pn.eic3 <- standardize.eic.fun(Pn.eic, eic.convergence.scale$rmtl)

    if (length(outcome.index)>1 && length(target)==length(outcome.index) &&
        !target.S && !is.null(eic.convergence.scale$rmst.free)) {
      Pn.eic3[[length(Pn.eic3)+1]] <-
        -Reduce("+", Pn.eic) /
        ((ifelse(any(eic.convergence.scale$rmst.free==0),
                 eic.convergence.scale$rmst.free+0.001,
                 eic.convergence.scale$rmst.free)) * sqrt(n))
      check.sup.norm <- max(abs(unlist(Pn.eic3)))<=(criterion)
    } else {
      check.sup.norm <- max(abs(unlist(Pn.eic3)))<=(criterion)
    }
    if (verbose) print(paste0("target level = ", criterion))
    if (verbose) print(paste0("max eic = ", max(abs(unlist(Pn.eic3)))))
    if (trace.eic.scales) {
      initial.scale <- list(rmtl=init.ic.rmtl, rmst.free=init.se.rmst.free)
      current.scale <- if (eic.scale=="current") eic.convergence.scale else eic.convergence.scale.fun(mat)
      print(paste0("max eic initial-scale = ", max.standardized.eic.fun(Pn.eic, initial.scale)))
      print(paste0("max eic current-scale = ", max.standardized.eic.fun(Pn.eic, current.scale)))
    }

    if (no.update.if.solved) {

      if (length(target)>1) {
        Pn.eic2 <- lapply(1:length(Pn.eic2), function(kk) {
          tmp <- Pn.eic2[[kk]]
          tmp[abs(Pn.eic3[[kk]])<=criterion] <- 0
          return(tmp)
        })
      } else {
        Pn.eic2[abs(Pn.eic3)<=criterion] <- 0
      }

      if (all(unlist(Pn.eic2)==0)) check.sup.norm <- TRUE else Pn.eic.norm <- Pn.eic.norm.fun(Pn.eic2, Pn.eic)
      print(Pn.eic2)

    }

    if (check.sup.norm | step==no.small.steps) {
      if (step==no.small.steps) {
        message("Warning: Algorithm did not converge")
      }

      if (verbose) print(paste0("converged", " at ", step, "th step"))
      if (verbose) print(paste0("eic = ", Pn.eic.fun(mat)))

      if (output.mat=="updated") {
        return(mat)
      }

      final.fit <- lapply(outcome.index[target], function(each) {
        sapply(1:length(tau), function(kk) {
          mean(rowSums(sapply(a, function(aa)
            (2*(aa==a[1])-1)*(mat[get(A.name)==aa,
                                  get(paste0("RMTL", estimation[[each]][["event"]],
                                             ".tau", kk))[1],
                                  by="id"][,2][[1]]))))
        })
      })

      names(final.fit) <- paste0("RMTL", sapply(outcome.index[target], function(each) estimation[[each]][["event"]]))
      final.ic <- eval.ic.rmtl(mat, final.fit, target.index=outcome.index[target])
      final.ic.rmtld.full <- eval.ic.rmtl(mat, final.fit,
                                          target.index=outcome.index[target],
                                          return.ic=TRUE)

      final.list <- lapply(1:length(final.fit), function(each.index) {
        out <- rbind(tmle.est=final.fit[[each.index]],
                     tmle.se=final.ic[[each.index]])
        colnames(out) <- paste0("tau=", tau)
        return(out)
      })

      names(final.list) <- paste0("RMTL", sapply(outcome.index[target], function(each) estimation[[each]]["event"]))

      tmle.list$tmle <- final.list
      if (simultaneous.ci) {
        final.sim.band.rmtld <- compute.simultaneous.rmtld.band(final.ic.rmtld.full,
                                                                final.fit,
                                                                target.index=outcome.index[target])
        tmle.list$q.max.95.rmtld <- data.table::rbindlist(lapply(outcome.index[target], function(each) {
          fit.delta <- estimation[[each]][["event"]]
          data.table::data.table(event=fit.delta,
                                 q.max.95=as.numeric(final.sim.band.rmtld[[paste0("RMTL", fit.delta)]]$q.max.95))
        }))
      } else {
        final.sim.band.rmtld <- NULL
      }
      tmle.list$tmle.rmtld <- format.tmle.rmtld.output(final.fit, final.ic,
                                                       sim.band=final.sim.band.rmtld,
                                                       target.index=outcome.index[target])
      final.fit.arm <- summarize.arm.rmtl.fit(mat, target.index=outcome.index[target])
      final.ic.arm.raw <- eval.ic.rmtl.arm(mat, final.fit.arm,
                                           target.index=outcome.index[target],
                                           return.ic=TRUE)
      final.ic.arm <- summarize.arm.se.from.ic(final.ic.arm.raw)
      final.fit.cif.arm <- summarize.arm.cif.fit(mat, target.index=outcome.index[target])
      final.ic.cif.arm <- eval.ic.cif.arm(mat, final.fit.cif.arm, target.index=outcome.index[target])
      tmle.list$tmle.group.rmtl <- format.arm.rmtl.output(final.fit.arm, final.ic.arm,
                                                          target.index=outcome.index[target])
      tmle.list$tmle.group.cif <- format.arm.cif.output(final.fit.cif.arm, final.ic.cif.arm,
                                                        target.index=outcome.index[target])
      tmle.list$tmle.adjusted.cif <- summarize.adjusted.cif(mat, target.index=outcome.index[target])
      if (length(outcome.index)>1 && length(target)==length(outcome.index) && !target.S) {
        final.fit.rmst.free <- list("RMST.free" = -Reduce("+", final.fit))
        final.ic.rmst.free.mat <- combine.rmst.free.ic(final.ic.rmtld.full)
        final.ic.rmst.free <- list("RMST.free" = sqrt(colMeans(final.ic.rmst.free.mat^2)/n))
        if (simultaneous.ci) {
          final.sim.band.rmst.free <- compute.simultaneous.scalar.band(
            final.ic.rmst.free.mat,
            final.fit.rmst.free[["RMST.free"]],
            seed.offset=3000
          )
          tmle.list$q.max.95.rmst.free <- data.table::data.table(
            q.max.95=as.numeric(final.sim.band.rmst.free$q.max.95)
          )
        } else {
          final.sim.band.rmst.free <- NULL
        }
        tmle.list$tmle.rmst.free <- format.tmle.rmst.free.output(
          final.fit.rmst.free,
          final.ic.rmst.free,
          sim.band=final.sim.band.rmst.free
        )

        final.fit.arm.rmst.free <- summarize.arm.rmst.free.fit(
          final.fit.arm,
          target.index=outcome.index[target]
        )
        final.ic.arm.rmst.free.raw <- combine.arm.rmst.free.ic(final.ic.arm.raw)
        final.ic.arm.rmst.free <- list(
          "RMST.free" = do.call("rbind", lapply(final.ic.arm.rmst.free.raw, function(ic.mat) {
            sqrt(colMeans(ic.mat^2)/n)
          }))
        )
        rownames(final.ic.arm.rmst.free[["RMST.free"]]) <- names(final.ic.arm.rmst.free.raw)
        colnames(final.ic.arm.rmst.free[["RMST.free"]]) <- paste0("tau=", tau)
        tmle.list$tmle.group.rmst.free <- format.arm.rmst.free.output(
          final.fit.arm.rmst.free,
          final.ic.arm.rmst.free
        )
      }

      if (check.sup.norm) {

        tmle.list$convergenced.at.step <- step

        break

      }
    }
  }

  return(tmle.list)

}
