prepare.data <- function(formula, data, subject, within){

if(!is.data.frame(data)){
  data <- as.data.frame(data)
}

dat <- model.frame(formula, data)
if (!(subject %in% names(data))){
  stop("The subject variable is not found!")
}
subject <- data[[subject]]
if (length(subject) != nrow(dat) && nrow(subject)!= nrow(data)){
  stop("There are missing values in the data.")
}

dat <- data.frame(dat, subject = subject)
nr_hypo <- attr(terms(formula), "factors")
perm_names <- t(attr(terms(formula), "factors")[-1, ])
fac_names <- colnames(nr_hypo)
fac_names_original <- fac_names
fac_names_simple <- colnames(perm_names)

if(is.null(fac_names_simple)){
  fac_names_simple <- within
}


if(!all(within %in% fac_names)){
  stop(paste0("The within-subjects factor ", 
              within[which(!(within %in% fac_names_simple))], 
              " is not part of the formula."))
}

outcome_names <- rownames(nr_hypo)[1]  # names of outcome variables
# extract names of outcome variables
if (grepl("cbind", outcome_names)){
  split1 <- strsplit(outcome_names, "(", fixed = TRUE)[[1]][-1]
  split2 <- strsplit(split1, ")", fixed = TRUE)[[1]]
  split3 <<- strsplit(split2, ",")[[1]]
} else {
  split3 <- outcome_names
}

EF <- rownames(nr_hypo)[-1]  # names of influencing factors
nf <- length(EF)
names(dat) <- c("response", EF, "subject")
#no. dimensions
p <- ncol(as.matrix(dat$response))
fl <- NA
for (aa in 1:nf) {
  fl[aa] <- nlevels(as.factor(dat[, (aa + 1)]))
}
levels <- list()
for (jj in 1:nf) {
  levels[[jj]] <- levels(as.factor(dat[, (jj + 1)]))
}
lev_names <- expand.grid(levels)
# number of hypotheses
tmp <- 0
for (i in 1:nf) {
  tmp <- c(tmp, choose(nf, i))
  nh <- sum(tmp)
}

names(fl) <- fac_names_simple

# determine within- and between-subject factors
no.subf <- length(within)
no.whole <- nf-no.subf
whole <- fac_names_simple[which(!fac_names_simple %in% within)]
lev.sub <- prod(fl[within])

# correct formula?
if (length(fac_names) != nf && length(fac_names) != nh){
  stop("Something is wrong with the formula. Please specify all or
       no interactions in crossed designs.")
}

# check that subjects are correctly labeled / no missing values / no within-factors ignored
if(nrow(data)/length(unique(subject)) != prod(fl[within])){
  stop(paste0("The number of subjects (", length(unique(subject)), ") times the 
  number of within-subject factor levels (", prod(fl[within]), ") does not equal
  the total number of observations (", nrow(data), "). 
  Check for missing values in the data."))
}


# check whether data is ordered correctly
if(no.whole >0 & !(whole[1] == colnames(dat)[2])){
  stop("The within-subjects factor(s) must be last in the formula!")
}


if (nf == 1) {
  # one-way layout
  dat2 <- dat[order(dat[, 2]), ]
  dat2 <- dat2[order(dat2[, "subject"]), ]
  #fac.groups <- dat2[, 2]
  hypo <- list(diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl))
  Y <- list(dat2)
  lev.sub <- fl
  names(fl) <- within
  # end one-way layout ------------------------------------------------------
} else {
  if(no.whole!=0){
    dat2 <- dat[do.call(order, dat[, c(2:(no.whole+1),      #order whole-plot factors
                                       ncol(dat),           # order by subject
                                       ((no.whole+2):(no.whole+2+no.subf-1)))]), ]# order sub-plot factors
  } else {
    dat2 <- dat[do.call(order, dat[, c(ncol(dat),           # order by subject
                                       (2:(2+no.subf-1)))]), ] # order sub-plot factors
  }
  lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
  if(length(whole) ==0){
    Y <- list(dat2)
  } else {
    Y<- split(dat2, dat2[, whole], lex.order = TRUE)
  }
}
nind <- sapply(Y, nrow)/lev.sub

## adapting formula argument, if interaction term missing
if (nrow(perm_names) != nh) {
  form2 <- as.formula(paste(outcome_names, "~", paste(fac_names, collapse = "*")))
  perm_names2 <- t(attr(terms(form2), "factors")[-1, ])
  fac_names2 <- attr(terms(form2), "term.labels")
  hyps <- HC(fl, perm_names2, fac_names2)
  hypo_matrices <- hyps[[1]]
  fac_names2 <- hyps[[2]]
  # choose only relevant entries of the hypo matrices
  indices <- grep(":", fac_names2, invert = TRUE)
  hypo <- lapply(indices, function(x) hypo_matrices[[x]])
  
} else if(nf !=1){
  hyp <- HC(fl, perm_names, fac_names)
  hypo <- hyp[[1]]
  fac_names <- hyp[[2]]
}

hypo_matrices <- lapply(hypo, function(x) x %x% diag(p))
# correcting for "empty" combinations (if no interaction specified)
n.groups <- prod(fl[whole])
if(nf != 1 & length(Y) != n.groups){
  index <- NULL
  for(i in 1:length(Y)){
    if(nrow(Y[[i]]) == 0){
      index <- c(index, i)
    }
  }
  Y <- Y[-index]
}

Ywhole <- lapply(Y, function(x) x$response)
if (p==1){
  Ywhole <- lapply(Ywhole, function(x) as.matrix(x))
}

Yw2 <- lapply(Ywhole, function(x) matrix(t(x), nrow = nrow(x)/prod(fl[within]), 
                                         ncol = p*prod(fl[within]), byrow=TRUE))

# ---------------------- error detection ------------------------------------
# no factor combinations with less than 2 observations
if (0 %in% nind || 1 %in% nind) {
  stop("There is at least one factor-level combination
           with less than 2 observations!")
}

out <- list(hypo = hypo_matrices, fac = fac_names, data = Yw2, n = nind, fl = fl,
           lev_names = lev_names, dat2 = dat2, levels = levels, no.subf = no.subf,
           EF = EF, nf = nf, fac_names_original = fac_names_original,
           no.whole = no.whole, lev.sub = lev.sub, p = p, split3 = split3)
return(out)

}
