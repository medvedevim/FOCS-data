rfile="data/D16 - 1000 nm_17.02.17_15.02.20.dat"

to.read = file(rfile, "rb")

# counting number of frames in file
nframes=file.info(rfile)$size/1024

# read file to array of unsigned integers
tdata<-readBin(to.read, integer(),size=2, n = nframes*512, endian = "little", signed = FALSE)

# reshape data to 512xN matrix
dim(tdata) <- c(512, nframes)
tdata <- t(tdata)

