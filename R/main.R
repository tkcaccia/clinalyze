


intersect = function (x, y, ...) 
{
  args <- list(x,y,...)
  n=length(args)
  inte=names(which(table(unlist(args))==n))
  if(length(inte)==0) inte=NULL
  inte
}

as.data.matrix = function(x){
  cn=colnames(x)
  rn=rownames(x)
  x=matrix(as.numeric(as.matrix(x)),ncol=ncol(x))
  colnames(x)=cn
  rownames(x)=rn
  x
}


continuous.test =
function (name, x, y, 
          center = c("median", "mean"),
          digits = 3, 
          scientific = FALSE, 
          range = c("IQR", "95%CI", "range", "sd"), 
          logchange = FALSE, pos = 2,
          method = c("non-parametric",    "parametric"), 
          total.column = FALSE, ...) 
{
    matchFUN.test = pmatch(method[1], c("non-parametric", "parametric"))
    if (matchFUN.test != 1 & matchFUN.test != 2) {
        stop("Method argument should one of \"non-parametric\",\"parametric\"")
    }
    y = as.factor(y)
    ll = levels(y)
    A = x[y == ll[1]]
    B = x[y == ll[2]]
    nn = length(levels(y))
    nn2 = length(unique(y[!is.na(x)]))
    v = data.frame(matrix(nrow = 1, ncol = nn + 3))
    v[1, 1] = name
    if (nn == 2 & nn2 > 1) {
        if (matchFUN.test == 1) {
            pval = wilcox.test(x ~ y, exact = FALSE, ...)$p.value
        }
        if (matchFUN.test == 2) {
            pval = t.test(x ~ y, ...)$p.value
        }
        if (logchange == TRUE) {
            fc = -log2(mean(A, na.rm = TRUE)/mean(B, na.rm = TRUE))
        }
    }
    if (nn > 2 & nn2 > 1) {
        if (matchFUN.test == 1) {
            pval = kruskal.test(x ~ y, ...)$p.value
        }
        if (matchFUN.test == 2) {
            pval = summary.aov(aov(x ~ y, ...))[[1]]$`Pr(>F)`[1]
        }
        logchange = FALSE
    }
    if (nn > 1) {
        v[1, 2:(1 + nn)] = tapply(x, y, function(x) txtsummary(x, f=center,
                                                               digits = digits, scientific = scientific, range = range))
        v[1, nn + 2] = txtsummary(x, f=center,digits = digits, scientific = scientific)
        if (nn2 == 1) {
            v[1, nn + 3] = NA
        }
        else {
            v[1, nn + 3] = format(pval, digits = 3, scientific = TRUE)
        }
    }

    if (pos == 1) {
            names(v) = c("Feature", paste(levels(y), ", ",center[1]," [",range[1],"]",sep = ""), 
                         paste("Total, ",center[1]," [",range[1],"]",sep = ""),
                               "p-value")  
    }
    else {
        
        v[1, 1] = paste(name,", ",center[1]," [",range[1],"]",sep = "")
                            
        names(v) = c("Feature", levels(y), "Total", "p-value")
    }
    v[v == "NA [NA NA]"] = "-"
    if (logchange == TRUE) {
        v = cbind(v[1, 1:(nn + 2)], logchange = round(fc, digits = 2), 
                  `p-value` = v[1, (nn + 3)])
        attr(v, "p-logchange") = fc
    }
    if (!total.column) {
        v = v[, -(nn + 2)]
    }
    if(length(x)<5000){                                
      sha=shapiro.test(x)$p.value 
      attr(v, "shapiro test") = sha                              
    }
    attr(v, "p-value") = pval
    return(v)
} 

                                  
multi_analysis =
  function (data, y, FUN = c("continuous.test", "correlation.test"), 
            ...) 
  {
    matchFUN = pmatch(FUN[1], c("continuous.test", "correlation.test"))
    if (is.na(matchFUN)) 
      stop("The function to be considered must be  \"continuous.test\" or \"correlation.test\".")
    if (matchFUN == 1) {
      FUN = continuous.test
    }
    if (matchFUN == 2) {
      FUN = correlation.test
    }
    da = NULL
    pval = NULL
    for (i in 1:ncol(data)) {
      sel.na = !is.na(data[, i])
      if (sum(sel.na) > 5) {
        temp = FUN(name = colnames(data)[i], x = data[sel.na,i], y = y[sel.na], ...)
        da = rbind(da, temp)
        pval[i] = attr(temp,"p-value")
      }
      else {
        if (matchFUN == 1) {
          da = rbind(da, c(colnames(data)[i], NA, NA,NA,NA))
        }
        if (matchFUN == 2) {
          da = rbind(da, c(colnames(data)[i], NA, NA))
        }
        
        pval[i] = NA
      }
    }
    FDR = p.adjust(pval, method = "fdr")
    FDR = format(FDR, digits = 3, scientific = TRUE)
    da = cbind(da, FDR)
    da
  }





#' Textual Summary of Numeric Data
#'
#' Produces a formatted string summarizing center and spread of a numeric variable.
#' @param x A numeric vector.
#' @param f Summary function: "mean" or "median".
#' @param digits Number of digits to round to.
#' @param scientific Logical; whether to use scientific notation.
#' @param range Type of dispersion: "IQR", "95%CI", "range", or "sd".
#' @return A string with center and variability summary.
#' @export
txtsummary =
function (x, f = c("median", "mean"), digits = 0, scientific = FALSE, 
          range = c("IQR", "95%CI", "range", "sd")) 
{
  matchFUN = pmatch(range[1], c("IQR", "95%CI", "range", "sd"))
  if (is.na(matchFUN)) 
    stop("The range to be considered must be \"IQR\", \"95%CI\", \"range\", or \"sd\".")
  matchf = pmatch(f[1], c("median", "mean"))
  if (is.na(matchFUN)) 
    stop("f to be considered must be \"median\" or \"mean\".")
  if (matchf == 1) 
    m = median(x,  na.rm = TRUE)
  if (matchf == 2) 
    m = mean(x, na.rm = TRUE)

  if (matchFUN == 1) 
    ci = quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  if (matchFUN == 2) 
    ci = quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  if (matchFUN == 3) 
    ci = range(x, na.rm = TRUE)
  if (matchFUN == 4) 
    ci = sd(x, na.rm = TRUE)
  if (scientific) {
    m = format(m, digits = digits, scientific = scientific)
    ci = format(ci, digits = digits, scientific = scientific)
  }
  else {
    m = round(m, digits = digits)
    ci = round(ci, digits = digits)
  }
  if (matchFUN == 4) {
    txt = paste(m, " [", ci, "]", sep = "")
  }
  else {
    txt = paste(m, " [", ci[1], " ", ci[2], "]", sep = "")
  }
  txt
}


                              
#-----for categorical data 
categorical.test =
function (name, x, y, total.column = FALSE, remove = "", ...) 
{
  if (!is.null(remove)) {
    x[x == remove] = NA
  }
  y = as.factor(y)
  nn = length(levels(y))
  t0 = table(x, y)
  ta = cbind(t0, as.matrix(table(x)))
  tb = sprintf("%.1f", t(t(ta)/colSums(ta)) * 100)
  tc = matrix(paste(ta, " (", tb, ")", sep = ""), ncol = nn + 
                1)
  tc[, c(colSums(t0), -1) == 0] = "-"
  v = NULL
  if (nrow(t0) == 1) {
    p.value = NA
    v[nn + 3] = ""
  }
  else {
    if(is.ordered(x) | is.ordered(y)){

      if(length(levels(y))>2 & length(levels(x))>2){
        p.value = cor.test(as.numeric(x),as.numeric(y),method="spearman")$p.value
      } else{ 
        p.value = jonckheere.test(as.numeric(x), as.numeric(y), alternative = "two.sided", nperm = 10000)$p.value
      }
    }else{
      p.value = fisher.test(t0, workspace = 10^7, ...)$p.value
    }
    v[nn + 3] = format(p.value, digits = 3, scientific = TRUE)
  }
  v[1] = name
  group = paste("   ", rownames(ta), ", n (%)", sep = "")
  cc = cbind(group, tc, rep(NA, length(group)))
  cc = rbind(v, cc)
  te = c(colnames(t0), "Total")
  colnames(cc) = c("Feature", te, "p-value")
  cc[is.na(cc)] = ""
  if (!total.column) {
    cc = cc[, -(nn + 2)]
  }
  attr(cc, "p-value") = p.value
  attr(cc, "shapiro test") = NA
  return(cc)
}






correlation.test= function(x,y,method = c("pearson", "spearman","MINE"), name=NA, perm=100 , ...){
  matchFUN = pmatch(method[1], c("pearson", "spearman","MINE"))
  if (is.na(matchFUN)) 
    stop("The method to be considered must be  \"pearson\", \"spearman\" or \"MINE\".")
  res=list()
  sel=!is.na(x) & !is.na(y)
  x=x[sel]
  y=y[sel]
  
  text = data.frame(matrix(nrow=1,ncol=3))
  text[1, 1] = name
  text[1,2]=NA
  text[1,3]=NA
  if(length(x)<5){
    warning("The number of correlated elements is less than 5.")
    estimate=NA
    p.value=NA
    
  }else{
    if(matchFUN==1){
      temp=cor.test(x,y,method="pearson")
      estimate=temp$estimate
      p.value=temp$p.value
    }
    if(matchFUN==2){
      temp=cor.test(x,y,method="spearman")
      estimate=temp$estimate
      p.value=temp$p.value
    }
    if(matchFUN==3){
      estimate=mine(x,y)$MIC
      v=NULL
      for(i in 1:perm){
        v[i]=mine(x,sample(y))$MIC
      }
      p.value=pnorm(estimate, mean = mean(v), sd = sqrt(((length(v) - 
                                                                    1)/length(v)) * var(v)), lower.tail = FALSE)
    }
    text[1,2]=round(estimate,digits=2)
    text[1,3]=format(p.value, digits = 3, scientific = TRUE)
  }
  if(matchFUN==1){
    names(text)=c("Feature","r","p-value")
  }  
  if(matchFUN==2){
    names(text)=c("Feature","rho","p-value")
  }  
  if(matchFUN==3){
    names(text)=c("Feature","MIC","p-value")
  }
  attr(text,"estimate")=estimate
  attr(text,"p-value")=p.value
  return(text)
}






frequency_matching =
  function (data, label, times = 5, seed = 1234) 
  {
    data = as.data.frame(data)
    data2 = data
    for (i in 1:ncol(data2)) {
      if (is.numeric(data2[, i])) {
        v <- quantile(data2[, i], prob = seq(0, 0.99, 0.2))
        data2[, i] = findInterval(data2[, i], v)
      }
      data2[, i] = as.vector(data2[, i])
      data2[is.na(data2[, i]), i] = "NA"
      data2[, i] = as.factor(data2[, i])
    }
    if (is.null(rownames(data2))) {
      rownames(data2) = paste("S", 1:nrow(data2), sep = "")
      rownames(data) = paste("S", 1:nrow(data), sep = "")
    }
    names(label) = rownames(data2)
    data2 = as.matrix(data2[!is.na(label), ])
    label = label[!is.na(label)]
    tal = table(label)
    minor = names(which.min(tal))
    labels_to_match = names(tal)[names(tal) != minor]
    data_minor = data2[label == minor, , drop = FALSE]
    nc = ncol(data2)
    selection = NULL
    for (lll in 1:length(labels_to_match)) {
      major = labels_to_match[lll]
      data_major = data2[label == major, , drop = FALSE]
      grid = list()
      count = list()
      rest = list()
      for (j in 1:nc) {
        lis = list()
        h = 1
        for (i in j:nc) {
          lis[[h]] = levels(as.factor(data2[, i]))
          h = h + 1
        }
        grid[[j]] = as.matrix(expand.grid(lis))
        co = apply(grid[[j]], 1, function(y) sum(apply(as.matrix(data_minor)[, 
                                                                             j:nc, drop = FALSE], 1, function(x) all(y == 
                                                                                                                       x))))
        count[[j]] = co * times[major]
      }
      rest = list()
      rest[[1]] = count[[1]]
      selected = rep(FALSE, nrow(data_major))
      names(selected) = rownames(data_major)
      for (j in 1:nc) {
        if (sum(rest[[j]]) > 0) {
          for (i in 1:nrow(grid[[j]])) {
            if (rest[[j]][i] != 0) {
              who = apply(as.matrix(data_major[, j:nc]), 
                          1, function(x) all(grid[[j]][i, ] == x))
              n_who = min(sum(who[!selected]), rest[[j]][i])
              rest[[j]][i] = rest[[j]][i] - n_who
              set.seed(seed)
              ss = sample(names(which(who[!selected])), 
                          n_who)
              selected[ss] = TRUE
            }
          }
          if (j < nc) {
            temp = list()
            for (ii in 2:ncol(grid[[j]])) temp[[ii - 1]] = as.matrix(grid[[j]])[, 
                                                                                ii]
            rest[[j + 1]] = aggregate(rest[[j]], by = temp, 
                                      FUN = sum, na.rm = TRUE)[, "x"]
          }
          else {
            rest[[j + 1]] = sum(rest[[j]])
          }
        }
      }
      if (sum(rest[[j]]) > 0) {
        set.seed(seed)
        ss = sample(which(!selected), rest[[j + 1]])
        selected[ss] = TRUE
      }
      selection = c(selection, rownames(data_major[selected, 
                                                   , drop = FALSE]))
    }
    selection = c(selection, rownames(data_minor))
    data = data[selection, ]
    data2 = data2[selection, ]
    label = label[selection]
    return(list(data = data, label = label, selection = selection))
  }




                        
multi_test = function(x, labels, ...){
  name_features=colnames(x)
  da=NULL
  pval=NULL
  shapiro=NULL
  for(i in 1:length(name_features)){
    if(is.numeric(x[,i])){
      temp=continuous.test(name_features[i],x[,i],labels, ...)
    } else{
      temp=categorical.test(name_features[i],x[,i],labels, ...)      
    }
    pval[i]=attr(temp,"p-value")
    shapiro[i]=attr(temp, "shapiro test")

     da=rbind(da,temp) 
  }

     
  
  attr(da,"p-value")=pval
  attr(da,"shapiro test")=shapiro
  return(da)
}






# Define unions for allowed types
setClassUnion("matrixOrNULL", c("data.frame","matrix", "NULL"))
setClassUnion("numericOrFactor", c("numeric", "factor"))

# Define the S4 class
setClass(
  "clinical.table",
  slots = c(x = "matrixOrNULL", 
            y = "numericOrFactor",
            total.column = "logical"),
  validity = function(object) {
    if (!is.null(object@x)) {
      if (!is.numeric(object@x)) {
        return("Slot 'x' must be a numeric matrix or NULL.")
      }
      if (nrow(object@x) != length(object@y)) {
        return("Number of rows in 'x' must match length of 'y'.")
      }
    }
    TRUE
  }
)

# Show only slot x
setMethod("show", "clinical.table", function(object) {
  print(object@x)
})


initialization = function(y, total.column = FALSE){
  ma <- new("clinical.table",
            x = NULL,
            y = y,
            total.column = total.column)
  ma
}

add_analysis = function(ma,name,x){
  if(is.numeric(x)){
    temp=continuous.test(name,x,ma@y,total.column = ma@total.column)
  }
  if(is.factor(x)){
    temp=categorical.test(name,x,ma@y,total.column = ma@total.column)
  }
  ma@x=rbind(ma@x,as.data.frame(temp))
  ma
}





