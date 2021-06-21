#####Thanks for Qiaozhi, for his talanted codes###
plot_grid_smooths <- function( df ,data_clean, plot_titles, ncol=2, xlab_in, ylab_in, ylim_in) { 
  out <- by (data = df, INDICES = df$variable, FUN = function(m) {
    
    droplevels(m$variable)
    tbl <- table(m$value)
    m$value <- factor(m$value, levels=levels(data_clean[,as.character(m$variable[1])]))
    plot_sub <- droplevels(m[m$value %in% names(tbl)[tbl >=5000],,drop=FALSE])
    ggplot( plot_sub, aes(y=Systolic_Average, x=Age, colour=(value))  )  + theme(legend.title=element_blank()) + ylab(ylab_in) +xlab (xlab_in) + geom_smooth() + coord_cartesian(y=ylim_in) +facet_wrap(~variable, scales = "free", labeller=as_labeller(plot_titles))  
    
  })
  do.call(ggarrange, c(out, ncol = ncol))
}


plot_grid_smooths_dbp <- function( df , plot_titles, ncol=2, xlab_in, ylab_in, ylim_in) { 
  out <- by (data = df, INDICES = df$variable, FUN = function(m) {
    droplevels(m$variable)
    tbl <- table(m$value)
    m$value <- factor(m$value, levels=levels(data_clean[,as.character(m$variable[1])]))
    plot_sub <- droplevels(m[m$value %in% names(tbl)[tbl >=5000],,drop=FALSE])
    ggplot( plot_sub, aes(y=Diastolic_Average, x=Age, colour=(value))  )  + theme(legend.title=element_blank()) + ylab(ylab_in) +xlab (xlab_in) + geom_smooth() + coord_cartesian(y=ylim_in) +facet_wrap(~variable, scales = "free", labeller=as_labeller(plot_titles))  
  })
  do.call(ggarrange, c(out, ncol = ncol))
}

plot_grid_smooths_bmi <- function( df , plot_titles, ncol=2, xlab_in, ylab_in, ylim_in) { 
  out <- by (data = df, INDICES = df$variable, FUN = function(m) {
    droplevels(m$variable)
    tbl <- table(m$value)
    m$value <- factor(m$value, levels=levels(data_clean[,as.character(m$variable[1])]))
    plot_sub <- droplevels(m[m$value %in% names(tbl)[tbl >=5000],,drop=FALSE])
    ggplot( plot_sub, aes(y=BMI, x=Age, colour=(value))  )  + theme(legend.title=element_blank()) + ylab(ylab_in) +xlab (xlab_in) + geom_smooth() + coord_cartesian(y=ylim_in) +facet_wrap(~variable, scales = "free", labeller=as_labeller(plot_titles))  
  })
  do.call(ggarrange, c(out, ncol = ncol))
}


plot_grid_smooths_pp <- function( df , plot_titles, ncol=2, xlab_in, ylab_in, ylim_in) { 
  out <- by (data = df, INDICES = df$variable, FUN = function(m) {
    droplevels(m$variable)
    tbl <- table(m$value)
    m$value <- factor(m$value, levels=levels(data_clean[,as.character(m$variable[1])]))
    plot_sub <- droplevels(m[m$value %in% names(tbl)[tbl >=5000],,drop=FALSE])
    ggplot( plot_sub, aes(y=pp, x=Age, colour=(value))  )  + theme(legend.title=element_blank()) + ylab(ylab_in) +xlab (xlab_in) + geom_smooth() + coord_cartesian(y=ylim_in) +facet_wrap(~variable, scales = "free", labeller=as_labeller(plot_titles))  
  })
  do.call(ggarrange, c(out, ncol = ncol))
}



spearman_confint <- function(s, n ) {
  
  lower <- tanh(atanh(s) -1.96/sqrt(n-3))
  upper <-tanh(atanh(s) +1.96/sqrt(n-3))
  confint <- upper - lower
  paste(round(s,2),sprintf('%.2f',round(confint,2)), sep="+/-")
}


coef_confint <- function(outlm, voi ) {
  coef <- coef(outlm)[voi]
  confint <- confint.default(outlm)[voi,2] - confint.default(outlm)[voi,1]
  paste(round(coef,2),round(confint,2), sep="+/-")
}

plot_grid_heats <- function( data_sub, text_var, ncol) { 
  df <- melt(data_sub[,c("TG", "HDL", text_var)], id =c("TG", "HDL"))
  out <- by (data = data_sub, INDICES = eval(parse(text=sprintf('data_sub$%s', text_var))), FUN = function(m) {
    ggplot( m, aes(x=TG, y=HDL))+ ggtitle(droplevels(eval(parse(text=sprintf('m$%s', text_var))))) +  stat_binhex() + theme(legend.position="none")+   scale_fill_gradient(low="#ffe5e5", high = "#ff0000") 
  })
  do.call(ggarrange, c(out, ncol = ncol))
}



# L2 dist -----------------------------------------------------------------

l2_dist <- function (x,y,z){
  
  
  z_ind <- which(z)
  z_nind <- which(!z)
  
  #z_ind <- sample(z_ind,size = min(length(z_ind), 1E3))
  #z_nind <- sample(z_nind,size = min(length(z_nind), 1E3))
  
  x[z_ind] <- scale(x[z_ind],scale = FALSE);
  y[z_ind] <- scale(y[z_ind],scale = FALSE);
  x[z_nind] <- scale(x[z_nind],scale = FALSE);
  y[z_nind] <- scale(y[z_nind],scale = FALSE);
  
  x_cut <- cut(x,breaks = 20)
  y_cut <- cut(y,breaks = 20)
  #x_cut <- cut(x,breaks = quantile(x,probs=seq(0,1,0.05)))
  #y_cut <- cut(y,breaks = quantile(y,probs = seq(0,1,0.05)))
  #y_cut <- cut(y,breaks = 20)
  count_1 <- table(x_cut[z_ind],y_cut[z_ind])
  dens_1 <- sweep(count_1,2, colSums(count_1),'/')
  dens_1[is.na(dens_1)] <- 0
  
  count_2 <- table(x_cut[z_nind],y_cut[z_nind])
  dens_2 <- sweep(count_2,2, colSums(count_2),'/')
  dens_2[is.na(dens_2)] <- 0
  
  val <- sqrt(sum((dens_1 - dens_2)^2))
}
l2_dist_test <- function(x,y,z,nits){
  val <- l2_dist(x,y,z)
  nulldist <-rep (0,nits);
  for (i in 1:nits){
    randvec <- rep(FALSE,length(x))
    randsamp <- sample.int(length(x), sum(z));
    randvec[randsamp] <- TRUE
    nulldist[i] <- l2_dist(x,y,randvec)
  }
  pval <- 1-pnorm(val, mean(nulldist), sd(nulldist))
  c(val,pval)
}



#
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}
# ----- END plot function ----- #