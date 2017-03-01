require(tools)
require(stringr)
require(ggplot2)
require(pspline)

# init variables
n = 25 # sliding average window

load_reflection_data <- function(datadir, n){
  
  # function to delete every N-th row
  Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
  
  # setup working directory and reading the files list
  setwd(datadir)
  
  tempdata <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  gdata <- data.frame(matrix(NA, nrow = 0, ncol = 4))
  
  l<-list.files(path=".",pattern = "*.dat$")
  
  for (rfile in l)
  {
    
    # get data about the sensor from the filename
    alloy <- str_extract(rfile, '.+?(?= - )')
    thickness <- str_extract(rfile, '[0-9]+?(?= nm)')
    
    ######### START OF BINARY READ ###############
    
    #read binary
    to.read = file(rfile, "rb")
    
    # counting number of frames in file
    nframes=file.info(rfile)$size/1024
    
    # read file to array of unsigned integers
    tdata<-readBin(to.read, integer(),size=2, n = nframes*512, endian = "little", signed = FALSE)
    
    # reshape data to 512xN matrix
    dim(tdata) <- c(512, nframes)
    tdata <- t(tdata)
    
    close(to.read)
    
    ######### END OF BINARY READ #################
    
    ######### START OF TEXT READ #################
    
    # read the according *.txt file by line and extract date-time using regexp
    con=file(paste(
      file_path_sans_ext(basename(rfile)),
      ".txt",
      sep=""),
      open="r")
    t<-str_extract_all(readLines(con), '\\d{1,2}.\\d{1,2}.\\d{2} \\d{1,2}:\\d{1,2}:\\d{1,2}', simplify = TRUE)
    
    #remove empty lines and convert to POSIX
    t <- t[!apply(t == "", 1, all),]
    t <- as.POSIXct(strptime(t, "%d.%m.%y %H:%M:%S"))
    
    #calculate time difference in seconds
    rtime <- as.numeric(difftime(t[2],t[1], units = "secs"))
    
    # create time sequence for the dataframe
    ttime<-seq(0,rtime,rtime/nframes)
    ttime<-as.matrix(ttime[1:nframes])
    
    close(con)
    
    ######### END OF TEXT READ #################
    
    # put together data and time
    tdata <- cbind(ttime, tdata)
    
    # assembling data frame
    tt<-(tdata[,100]-tdata[nframes,100])/(tdata[1,100]-tdata[nframes,100])
    tempdata<-cbind(tdata[,1],tt)
    tempdata <- aggregate(tempdata,list(rep(1:(nrow(tempdata)%/%n+1),each=n,len=nrow(tempdata))),mean)[-1];
    tempdata<-cbind(tempdata, rep(alloy, nrow(tempdata)), rep(thickness, nrow(tempdata)))
    tempdata<-as.data.frame(tempdata)
    gdata <- rbind(gdata,tempdata)
    
  }
  
  # preparing dataframe for publishing
  names(gdata)<-c("Time", "Reflection", "Alloy", "Thickness")
  gdata$Time <- as.numeric(as.character(gdata$Time))
  gdata$Reflection <- as.numeric(as.character(gdata$Reflection))
  gdata<-gdata[which(gdata$Reflection<1.1 & gdata$Reflection>=0),]
  gdata$Thickness <- factor(gdata$Thickness, levels = as.character(sort(as.numeric(levels(gdata$Thickness)))))
  
  # clean up
  rm(tempdata, l, rfile)
  return(gdata)
}
find_equilibrium <- function(dataframe, roundv){
  deriv <- predict(sm.spline(dataframe$Time, dataframe$Reflection), dataframe$Time, 1)
  max_deriv <- which.max(abs(deriv))
  tt <- dataframe[max_deriv:nrow(deriv),]
  return(round(tt[which(tt$Reflection<1e-03)[1],"Time"], roundv))
}

gdata <- load_reflection_data("C:/R/FOCS/data", n)

tmp<-gdata[which(gdata$Alloy=="D16"),]
tmp2<-tmp[which(tmp$Reflection<0.95),]
tmp3<-tmp[which(tmp$Reflection>=0.95),]

####### Plot for different alloys
for (i in levels(gdata$Alloy))
{
  tmp <- gdata[which(gdata$Alloy==i),]
  print(ggplot(aes(x=Time, y=Reflection, shape=Thickness), data = tmp)+
    geom_point()+geom_path()+theme_bw()+
    theme(text = element_text(size=14))+ ggtitle(i))
}


eqtable <- data.frame(matrix(NA, nrow = 1, ncol = (length(levels(gdata$Thickness))+1)))
names(eqtable) <- c("Alloy", levels(gdata$Thickness))
for (i in levels(gdata$Alloy))
{
  foo <- c()
  for (f in levels(gdata$Thickness))
  {
    tt <- gdata[which(gdata$Thickness==f & gdata$Alloy==i),]
    foo <- cbind(foo, find_equilibrium(tt,0))
  }
  eqtable <- rbind(eqtable, c(i,as.numeric(foo)))
  names(eqtable) <- c("Alloy", levels(gdata$Thickness))
}

