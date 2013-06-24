polarx <- function(radius, theta) { radius * cos(theta) }
polary <- function(radius, theta) { radius * sin(theta) }

target <- function(moa, count) {
       r <- rnorm(count*100, 0, moa/2.0);
       r <- sample(r, count);
       a <- (sample(0:100000, count)/100000)*pi;
       p <- list(x=polarx(r,a), y=polary(r,a));
       m <- max(abs(c(p$x,p$y,1.5)));
       l <- c(-m,m);
       t <- c(0.5,1,2,3);
       plot.new();
       frame();
       plot.window(xlim=l,ylim=l);
       points(p,pch=16,col='#0000ff33');
       symbols(rep(0,length(t)),rep(0,length(t)),circles=t/2.0,add=T,lwd=5,inches=F,fg="#00000055")
       axis(1);
       title(main=sprintf("Target of %d shots with SD at %.2f-MOA", count, moa), xlab="Dispersion in MOA");
       list(x=p$x,y=p$y)
}

group <- function(l) {
      x <- l$x;
      y <- l$y;
      md <- 0;
      for(i in 1:length(x)) {
            for(j in 1:length(x)) {
                  dist <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2);
                  if(dist > md) md = dist;
            }
      }
      md
}

groups <- function(count,points) {
       sizes <- c();
       for(i in 1:(length(points$x) - count + 1)) {
             p <- list(x=points$x[seq(i,i+count-1)],
                       y=points$y[seq(i,i+count-1)]);
             sizes <- c(sizes,group(p));
       }
       sizes
}

testRifle <- function(moa, count=1000) {
  r <- target(moa,count)
  for(i in 2:20) { 
    gs <- groups(i,r);
    g <- mean(gs);
    s <- sd(gs);
    print(sprintf("%4g %4d shot groups, size: %.2f moa, stdev: %.2f moa, factor: %.2f", floor(count/i), i, g, s, moa/g));
  }
}
testRifle(0.7)

rsurf <- function(moa=1) {
    gran <- 100
    ncol <- 32
    totl <- 50000
    r <- target(moa,totl)
    #r$x <- r$x[r$x > -1.2*moa & r$x < 1.2*moa]
    #r$y <- r$y[r$y > -1.2*moa & r$y < 1.2*moa]
    sx <- split(r$x,cut(r$x, breaks=gran))
    sy <- split(r$y,cut(r$y, breaks=gran))
    pz <- 1:(gran*gran)
    count <- 1
    for(i in 1:gran) {
          for(j in 1:gran) {
                xs <- r$x %in% sx[[j]]
                ys <- r$y %in% sy[[i]]
                pz[count] <- length((1:gran)[xs & ys])
                count <- count + 1
          }
    }
    m <- matrix(pz, gran, gran)
    nr <- nrow(m)
    nc <- ncol(m)
    mf <- m[-1,-1] + m[-1,-nc] + m[-nr,-1] + m[-nr, -nc]
    cols <- cut(mf, ncol)
    ramp <- topo.colors(ncol)[cols]
    persp(seq(min(r$x), max(r$x), length=(gran)),
          seq(min(r$y), max(r$y), length=(gran)),
          m, theta = 20, phi = 20, box=T,
          xlab="Windage (MOA)", ylab="Elevation (MOA)",
          zlab="Density", expand=0.5,
          col=ramp, border="#00000033", bg="black")
    title(main=sprintf("%g shots at SD = %.2f MOA", totl, moa))
    cols
}
c <- rsurf()
