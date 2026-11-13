# Ruby examples

## Jabber bot

```
require 'xmpp4r/client'

class JabberBot < Jabber::Client
  include Jabber

  def initialize(jabber_id, jabber_password, options = {})
    @jabber_id = jabber_id
    @jabber_password = jabber_password
    @options = options
    @jid = JID::new(@jabber_id)
    # super calls the base_method
    super(@jid)
  end

  def on_message(from, body)
    raise "on_message not implemented!"
  end

  def connect_and_authenticate
    connect
    auth(@jabber_password)
    send(Presence::new)
    puts "Connected! Send messages to #{@jid.strip.to_s}."
  end

  def run
    main = Thread.current
    add_message_callback { |message|
      next if message.type == :error or message.body.blank?
      # support exit command if owner is set
      main.wakeup if @options[:owner_jabber_id] == message.from && message.body == 'exit'
      # simple callback
      self.on_message(message.from, message.body)
    }
    Thread.stop
    close
  end

  def say(jabber_id, body)
    send(Message::new(jabber_id, body).set_type(:chat).set_id('1'))
  end

end
```

## Codebits rfid bot

```
#!/usr/bin/env ruby
require 'jabber/bot'
require 'socket'
require 'date'

def sendreceive(request,ip='85.247.250.5',port=13568)
    puts "initializing socket..."
    TCPSocket.open(ip, port.to_i) do |sock|
        puts "done\n"
        puts "sending request..."
        sock.write(request)
        puts "done\n"
        puts "getting answer..."
        result = sock.read(100)
        puts "done\n"
        sock.close
        return result
    end
end

def timefromnow(string)
    now = DateTime.now
    date = DateTime.parse(string)
    days = (now-date).to_f
    secs = days * 86400
    if secs < 60
        msg = secs.to_i + " seconds ago."
    else
        if secs > 59 && secs < 120
            msg = "1 minute ago."
        else
            mins = secs / 60
            if mins < 60
                msg = mins.to_i.to_s + " minutes ago."
            else
                if mins > 59 && mins < 120
                    msg = "1 hour ago."
                else
                    hours = mins / 60
                    if hours < 24
                        msg = hours.to_i.to_s + " hours ago."
                    else
                        if hours > 23 && hours < 48
                            msg = "1 day ago."
                        else
                            days = hours / 24
                            msg = days.to_i.to_s + " days ago."
                        end
                    end
                end
            end
        end
    end

    return msg

end

class Rfidcodebits

    def initialize(file = "rfidbits.txt")
        @cache = file
        @results = Array.new
    end

    def whereis(message)

        flag = false
        File.open(@cache,'r') { |fs|
            pat=message.gsub(/ /,' .*')
            while line = fs.gets
                words = line.split(/\|/)
                if words[5].match(/#{pat}/i) || words[2].match(/^#{message.strip}$/)
                    flag = true
                    @results.push("#{words[5]} (#{words[2]}) was last seen at " + words[1] + ", " + timefromnow(words[3]))
                end
            end
            fs.close
            if !flag
                @results.push("#{message} wasn't found.")
            end
            @results
        }

    end

end

# Create a public Jabber::Bot
bot = Jabber::Bot.new(
  :jabber_id => 'codebits@jabber.cc',
  :password => 'hackathon',
  :master    => 'guillaume.riflet@gmail.com',
  :is_public => true
)

# Give your bot a public command
bot.add_command(
  :syntax      => 'rand',
  :description => 'Produce a random number from 0 to 10',
  :regex       => /^rand$/,
  :is_public   => false
) { rand(10).to_s }

# Give your bot a private command with an alias
bot.add_command(
  :syntax      => 'puts ',
  :description => 'Write something to $stdout',
  :regex       => /^puts\s+.+$/,
  :alias       => [
      :syntax => 'p ',
      :regex => /^p\s+.+$/
  ]
) do |sender, message|
  puts message
  "'#{message}' written to $stdout"
end

bot.add_command(
    :syntax         => 'whereis ',
    :description    => 'People spotter, give name or id',
    :regex          => /^whereis */i,
    :alias          => [
        :syntax => 'W',
        :regex  => /^W */i
    ],
    :is_public      => true
) do |sender, message|
    #sendreceive("whereis:#{message}")
    bits = Rfidcodebits.new
    bits.whereis(message).join("\n\n")
end

bot.add_command(
   :syntax      => 'kill',
   :description => 'Disconnect me, master.',
   :regex       => /^kill$/,
   :alias       => [
     :syntax => 'K',
     :regex => /^K$/
   ]
) {
  "Ah, you killed me master."
  self.disconnect
}

# Bring your new bot to life
begin
  bot.connect
rescue
  puts "The bot got, somehow, disconnected."
end

#puts sendreceive("Guillaume\n", 'localhost', 3000)
#puts sendreceive("Guillaume")

#mybits = Rfidcodebits.new
#puts mybits.whereis(" 32 ").join("\n")
#puts timefromnow('2008-11-14 17:18:46')
```

## PatternTransform

In this example, we'll see how to use a couple of classes in a ruby program in order to execute in a batch a series of actions to perform with the MOHID utility tools, such as [ConvertToHdf5](../../tools/converttohdf5.md) or Hdf5Extractor.

Concretely speaking, the example takes a template file with the configuration of the MOHID utility (such as *ConvertToHdf5Action.dat* or *Hdf5Extractor.dat*) and file with a list of substitutions to perform on the template configuration file, one line per configuration file.

### Template file

Here is an example of the template file:

```
#root_Hydro_HDF5Extractor.dat
OUTPUTFILENAME           : MOHID_HYDRO_XX.hdf5

FILENAME                 : ../Glue/MOHID_HYDRO_XX.hdf5

START_TIME               : 2004 MM SDD 0 0 0
END_TIME                 : 2004 MM EDD 0 0 0

CONVERT_V3_TO_V4         : 0

<BeginParameter>
PROPERTY                 : velocity U
HDF_GROUP                : /Results/velocity U
<EndParameter>

<BeginParameter>
PROPERTY                 : velocity V
HDF_GROUP                : /Results/velocity V
<EndParameter>

<BeginParameter>
PROPERTY                 : water level
HDF_GROUP                : /Results/water level
<EndParameter>
```

### Substitutions list

Here is the file containing the list of substitutions to perform on the above template:

```
#List.txt
(XX) (MM) (EDD) (SDD)
01    1     31  2
02    2     29  1
03    3     31  1
04    4     30  1
05    5     31  1
06    6     30  1
07    7     31  1
08    8     31  1
09    9     30  1
10    10    31  1
11    11    30  1
12    12    31  1
```

The first line contains the list of strings to match in the template file (such as *XX* or *MM*). Each string to match is between parentheses. The strings to match are space-separated. The following lines containg the strings to substitute, space-separated. Each new line will correspond to a new configuration file.

For example, the substitution of *(XX) (MM) (EDD) (SDD)* with *01 1 31 2* on the template file will yield the following result file:

```
#Hydro_HDF5Extractor_001.dat
OUTPUTFILENAME           : MOHID_HYDRO_01.hdf5

FILENAME                 : ../Glue/MOHID_HYDRO_01.hdf5

START_TIME               : 2004 1 2 0 0 0
END_TIME                 : 2004 1 31 0 0 0

CONVERT_V3_TO_V4         : 0

<BeginParameter>
PROPERTY                 : velocity U
HDF_GROUP                : /Results/velocity U
<EndParameter>

<BeginParameter>
PROPERTY                 : velocity V
HDF_GROUP                : /Results/velocity V
<EndParameter>

<BeginParameter>
PROPERTY                 : water level
HDF_GROUP                : /Results/water level
<EndParameter>
```

### Program file

Here is the program code that will generate the suite of configuration files and run the batch of transformations with *HDF5Extract*:

```
require "patterntransform"

hextract = IsTheActionToPerformWithFiles.new("root_Hydro_HDF5Extractor.dat", "list.txt")
hextract.transformsInto("Hydro_HDF5Extractor")
hextract.usesAGenericMohidtoolBatchfile("HydroExtract.bat", "HDF5Extractor.exe", "HDF5Extractor.dat")
hextract.ing

wextract = IsTheActionToPerformWithFiles.new("root_Water_HDF5Extractor.dat", "list.txt")
wextract.transformsInto("Water_HDF5Extractor")
wextract.usesAGenericMohidtoolBatchfile("WaterExtract.bat", "HDF5Extractor.exe", "HDF5Extractor.dat")
wextract.ing
```

### Class IsTheActionToPerformWith

```
class IsTheActionToPerformWith

  attr_reader :patterncontent, :listofsubs
  attr_reader :matchline, :batch

  def initialize (patterncontent, listofsubs)
    @patterncontent, @listofsubs = patterncontent, listofsubs
    @listofsubs = @listofsubs.to_a if !@listofsubs.is_a?(Array)
    @matchline = @listofsubs.shift
    @matches = @matchline.tr('()','').split(' ')
    @closings = @matchline.gsub(/\(.+?\)/, '*').split(/\*| /)
    puts @closings.join(' ')
  end

  #private
  def thePatternIsTransformedByTheBlock
    @j = 1
    @listofsubs.each do |line|
      yield line
      @j = @j + 1
    end
  end

  #private
  def buildHashOfSubs(array)
    hash = Hash.new
    array.size.times { |i| hash.store(@matches[i], @closings[2*i].chomp + array[i] + @closings[2*i + 1].chomp) }
    return hash
  end

  def transformsInto(filename)

    @outfilename = filename
    thePatternIsTransformedByTheBlock{ |line|
      filecontent = String.new( @patterncontent )

      tokens = line.split( ' ' )
      subhash = buildHashOfSubs(tokens)

      @matches.each { |element|
        iselement = false
        subhash.each{ |match,sub|
            if (match.eql? element)
                iselement = true
            end
        }
        if !iselement
            filecontent = String.new( filecontent.reject{|n| n.match(element)}.to_s )
        end
      }
      subhash.each { |match, sub| filecontent.replace filecontent.gsub(match, sub) }

      tokens = line.split( ' ' )
      subhash = buildHashOfSubs(tokens)
      subhash.each { |match, sub| filecontent.replace filecontent.gsub(match, sub) }

      mytrgfile = File.new(filename + "_" + @j.to_s.rjust(3,"0")  + ".dat", "w")
      mytrgfile.puts filecontent
      mytrgfile.close

      puts filecontent
    }

  end

  #Running
  def ing
    puts "Running #{@batchfilename}..."
  end

  #Running
  def ning
    ing
  end

#Method 1
  def usesTheFollowingBatchfileNameAndBlock(batchfilename)

    @batchfilename = batchfilename

    @batch = String.new

    batchhead = <<HEAD
    @echo off
    for %%i in (#{@outfilename}_*.dat) do (
HEAD

    yield

    batchfoot = <<FOOT
    )
    pause
    @echo on
FOOT

    batchfile = File.new(@batchfilename, "w")
    batchfile.puts batchhead + @batch + batchfoot
    batchfile.close

  end

  #Method 2
  def usesAGenericMohidtoolBatchfile(mohidtool, batchfilename, configfilename)
    usesTheFollowingBatchfileNameAndBlock(batchfilename){
      @batch = <<BATCH
      copy %%i #{configfilename}
      #{mohidtool}
BATCH
    }
  end

  def usesAGenericMohidtool(mohidtool, configfilename)
    usesAGenericMohidtoolBatchfile(mohidtool, "All_UseTool.bat", configfilename)
  end

#Methods 3
  def usesAConverttohdf5BatchfileNamed(batchfilename)
    usesAGenericMohidtoolBatchfile( "ConvertToHdf5.exe", batchfilename, "ConvertToHDF5Action.dat")
  end

  def usesAConvert2netcdfBatchfileNamed(batchfilename)
    usesAGenericMohidtoolBatchfile( "Convert2Netcdf.exe", batchfilename, "Convert2Netcdf.dat")
  end

end

class IsTheActionToPerformWithFiles < IsTheActionToPerformWith

  def initialize(patternfilename, listofsubsfilename)

    if File.exists?(patternfilename) and File.exists?(listofsubsfilename) then

      #patternfile = File.open(patternfilename)
      patterncontent = String.new( File.open(patternfilename).read )

      #listofsubsfile = File.open(listofsubsfilename)
      listofsubs = String.new( File.open(listofsubsfilename).read )

      super(patterncontent, listofsubs)

    else

      puts "Error: one or both of #{patternfilename} or #{listofsubsfilename} doesn't exist."

    end

  end

end

if __FILE__ == $0

mycontent = <<CONVERTTOHDF5ACTION
<begin_file>
ACTION                   : GLUES HDF5 FILES

OUTPUTFILENAME           : MOHID_WATER_01.hdf5

<<begin_listofsubs>>
..\\Extraction\\Stride_WaterProperties_1.hdf5
..\\Extraction\\Stride_WaterProperties_2.hdf5
<<end_listofsubs>>
<end_file>
CONVERTTOHDF5ACTION

mylistofsubs = <<listofsubs
_(01) s_(1). s_(2). (WATER) (WaterProperties)
   01      1        2    WATER   WaterProperties
   01      1        2    HYDRO   Hydrodynamic
   02      3        5    WATER   WaterProperties
   02      3        5    HYDRO   Hydrodynamic
listofsubs

  ## 1 Create ##########
  #You can create the class object with files ...
  glue = IsTheActionToPerformWithFiles.new("root_ConvertToHDF5Action.dat", "listofsubs.txt")
  #... or with strings
  glue = IsTheActionToPerformWith.new(mycontent, mylistofsubs)

  ## 2 Generate new files ##########
  glue.transformsInto("GlueAction")

  ## 3 Generate runscript ##########
  #You can easily customize the batchfile main loop ...
  glue.usesTheFollowingBatchfileNameAndBlock("GlueAllHdf5.bat") {
    glue.batch.replace <<BATCH
      copy %%i ConvertToHdf5Action.dat
      ConvertToHdf5.exe
BATCH
  }
  #... or use a standard one
  glue.usesAConverttohdf5BatchfileNamed("GlueAllHdf5.bat")

  ## 4 Run ##########
  glue.ing

end
```
