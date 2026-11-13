# MovAverageRuby

```
#This code computes the moving average from mohid type .xyz files.
def getData(arg1, arg2)
  if File.exist? arg1
    data = []
    File.open(arg1, "r") { |fs|
      while line = fs.gets
        tokens = line.split(" ")
	    if tokens.size == 3
  		  data.push(tokens[arg2])
	    end
      end
      fs.close
    }
    return data
  end
end

def makeMovAverage(time,interval,data)
  out = []
  n = time.size - 1
  time.each { |t|
    sum = 0
    count = 0
    for i in 0..n do
      if ( t.to_f - time[i].to_f ).abs < interval.to_f then
	    sum = sum.to_f + data[i].to_f
	    count = count.to_i + 1
	  end
    end
    out.push(sum/count)
  }
  return out
end

def writeMovAverage(lon,lat,data,filename)
  n = lon.size-1
  content = ["<begin_xyz>\n"]
  for i in 0..n do
    line = [lon[i],lat[i],data[i]]
    content.push( line.join(" ") + "\n" )
  end
  content.push("<end_xyz>\n")
  File.open(filename, "w") { |fw|
    fw.puts content
    fw.close
  }
  puts "Done #{filename}"
end

def computeMovAverage(lon,lat,time,interval,filename)
  data =  getData(filename, 2)
  newData = makeMovAverage(time, interval,data)
  writeMovAverage(lon,lat,newData,filename.split('_').first + "_MOVAV.xyz")
end

lon = getData("table1Sec_AB.xyz", 0)
lat = getData("table1Sec_AB.xyz", 1)
time = getData("table1Sec_AB.xyz", 2)
puts "Time: #{time.size}"

computeMovAverage(lon,lat,time,120,"table1WindEast_AB.xyz")
computeMovAverage(lon,lat,time,120,"table1WindNorth_AB.xyz")
computeMovAverage(lon,lat,time,120,"table1WindUp_AB.xyz")
```
