require(tools)
require(stringr)

# setup working directory and reading the files list
setwd(paste(getwd(),"/data", sep=""))

l<-list.files(path=".",pattern = "*.dat$")

rfile="D16 - 300 nm_17.02.17_14.55.37.dat"
to.read = file(rfile, "rb")



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

# counting number of frames in file
nframes=file.info(rfile)$size/1024

# create time sequence
ttime<-seq(0,rtime,rtime/nframes)
ttime<-as.matrix(ttime[1:nframes])

# read file to array of unsigned integers
tdata<-readBin(to.read, integer(),size=2, n = nframes*512, endian = "little", signed = FALSE)

# reshape data to 512xN matrix
dim(tdata) <- c(512, nframes)
tdata <- t(tdata)

tdata <- cbind(ttime, tdata)
