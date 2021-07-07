using Printf
using Plots
using StaticArrays

mutable struct ChannelStruct
    native_channel_name::String
    custom_channel_name::String
    native_order::Int16
    custom_order::Int16
    board_stream::Int16
    chip_channel::Int16
    port_name::String
    port_prefix::String
    port_number::Int16
    electrode_impedance_magnitude::Float32
    electrode_impedance_phase::Float32
end

mutable struct SpikeTriggerStruct
    voltage_trigger_mode::Int16
    voltage_threshold::Int16
    digital_trigger_channel::Int16
    digital_edge_polarity::Int16
end

mutable struct Version
    major::Int
    minor::Int
    Version() = new((0 for _ in 1:length(fieldnames(Version)))...)
    #Version(args...) = new(args...)
end

mutable struct Indices
    amplifier::Int
    aux_input::Int
    supply_voltage::Int
    board_adc::Int
    board_dig_in::Int
    board_dig_out::Int
    Indices() = new((1 for _ in 1:length(fieldnames(Indices)))...)
end

mutable struct Freq
    dsp_enabled::Int16
    actual_dsp_cutoff_frequency::Float32
    actual_lower_bandwidth::Float32
    actual_upper_bandwidth::Float32
    desired_dsp_cutoff_frequency::Float32
    desired_lower_bandwidth::Float32
    desired_upper_bandwidth::Float32
    notch_filter_frequency::Float32
    desired_impedance_test_frequency::Float32
    actual_impedance_test_frequency::Float32
    amplifier_sample_rate::Float32
    aux_input_sample_rate::Float32
    supply_voltage_sample_rate::Float32
    board_adc_sample_rate::Float32
    board_dig_in_sample_rate::Float32
    Freq() = new((0 for _ in 1:length(fieldnames(Freq)))...)
end

mutable struct Header
    version::Version
    frequency_parameters::Freq
    sample_rate::Float32
    num_temp_sensor_channels::Int16
    eval_board_mode::Int16
    reference_channel::String
    num_samples_per_data_block::Int16
    num_amplifier_channels::Int16
    num_aux_input_channels::Int16
    num_supply_voltage_channels::Int16
    num_board_adc_channels::Int16
    num_board_dig_in_channels::Int16
    num_board_dig_out_channels::Int16
    notes
    spike_triggers
    amplifier_channels
    aux_input_channels
    supply_voltage_channels
    temp_sensor_channels
    board_adc_channels
    board_dig_in_channels
    board_dig_out_channels
    Header() = new(Version(), Freq(), 0, 0, 0, "", 0, 0, 0, 0, 0, 0, 0)
end 

mutable struct Data
    t_amplifier
    t_aux_input
    t_supply_voltage
    t_board_adc
    t_dig
    t_temp_sensor
    amplifier_data
    aux_input_data
    supply_voltage_data
    temp_sensor_data
    board_adc_data
    board_dig_in_data
    board_dig_in_raw
    board_dig_out_data
    board_dig_out_raw
    Data() = new()
end

mutable struct Result
    t_amplifier
    t_aux_input
    t_supply_voltage
    t_board_adc
    t_dig
    t_temp_sensor
    spike_triggers
    notes
    frequency_parameters
    reference_channel
    amplifier_channels
    amplifier_data
    aux_input_channels
    aux_input_data
    supply_voltage_channels
    supply_voltage_data
    board_adc_channels
    board_adc_data
    board_dig_in_channels
    board_dig_in_data
    board_dig_out_channels
    board_dig_out_data
    Result() = new()
end


function readQString(fid)
    #= Read Qt style String. The first 32-bit unsigned number indicates the length of the string (in bytes).
    If this number equals 0xffffffff, the string is null =#
    
    a = ""
    length = read(fid, UInt32)
    if length == 0xffffffff
        return
    end
    
    # Convert length from bytes to 16-bit Unicode words
    length = length / 2
    for i = 1:length
        thisChar = Char(read(fid, UInt16))
        a = a * thisChar
    end
    return a
end


function plural(n)
    # s = plural(n)
    # Utility function to optionally pluralize words based on the value of n
    if n == 1
        s = ""
    else
        s = "s"
    end
    return s
end


function get_bytes_per_data_block(header)
    # Calculates the number of bytes in each 60 or 128 sample datablock
    
    # Each data block contains 60 or 128 amplifier samples
    bytes_per_block = header.num_samples_per_data_block * 4 # timestamp data
    bytes_per_block = bytes_per_block + header.num_samples_per_data_block * 2 * header.num_amplifier_channels
    
    # Auxiliary inputs are sampled 4x slower than amplifiers
    bytes_per_block = bytes_per_block + (header.num_samples_per_data_block / 4) * 2 * header.num_aux_input_channels
    
    # Supply voltage is sampled 60 or 128x slower than amplifiers
    bytes_per_block = bytes_per_block + 1 * 2 * header.num_supply_voltage_channels
    
    # Board analog inputs are sampled at same rate as amplifiers
    bytes_per_block = bytes_per_block + header.num_samples_per_data_block * 2 * header.num_board_adc_channels
    
    # Board digital inputs are sampled at same rate as amplifiers
    if header.num_board_dig_in_channels > 0
        bytes_per_block = bytes_per_block + header.num_samples_per_data_block * 2
    end
        
    # Board digital outputs are sampled at same rate as amplifiers
    if header.num_board_dig_out_channels > 0
        bytes_per_block = bytes_per_block + header.num_samples_per_data_block * 2
    end
    
    # Temp sensor is sampled 60 or 128x slower than amplifiers
    if header.num_temp_sensor_channels > 0
        bytes_per_block = bytes_per_block + 1 * 2 * header.num_temp_sensor_channels
    end
    
    return bytes_per_block
end

# Define read_header function
function read_header(fid)
    # Check 'magic number' at beginning of file to make sure this is an Intan Technologies RHD2000 data file.
    magic_number = read(fid, UInt32)
    if magic_number != 0xc6912702
        error("Unrecognized file type.")
    end
    
    header = Header()
    # Read version number
    version = Version()
    version.major = read(fid, Int16)
    version.minor = read(fid, Int16)
    
    println("\nReading Intan Technologies RHD2000 Data File, Version ", version.major, ".", version.minor)
    
    # Read information of sampling rate and amplifier frequency settings
    header.sample_rate = read(fid, Float32)
    header.frequency_parameters.dsp_enabled = read(fid, Int16)
    header.frequency_parameters.actual_dsp_cutoff_frequency = read(fid, Float32)
    header.frequency_parameters.actual_lower_bandwidth = read(fid, Float32)
    header.frequency_parameters.actual_upper_bandwidth = read(fid, Float32)
    header.frequency_parameters.desired_dsp_cutoff_frequency = read(fid, Float32)
    header.frequency_parameters.desired_lower_bandwidth = read(fid, Float32)
    header.frequency_parameters.desired_upper_bandwidth = read(fid, Float32)
    
    # This tells us if a software 50/60 Hz notch filter was enabled during the data acquisition.
    notch_filter_mode = read(fid, Int16)
    notch_filter_frequency = 0
    if notch_filter_mode == 1
        notch_filter_frequency = 50
    elseif notch_filter_mode == 2
        notch_filter_frequency = 60
    end
    header.frequency_parameters.notch_filter_frequency = notch_filter_frequency
    
    header.frequency_parameters.desired_impedance_test_frequency = read(fid, Float32)
    header.frequency_parameters.actual_impedance_test_frequency = read(fid, Float32)
    
    # Place notes in array of Strings
    header.notes = [readQString(fid), readQString(fid), readQString(fid)]
    
    # If data file is from GUI v1.1 or later, see if temperature sensor data was saved
    num_temp_sensor_channels = 0
    if (version.major == 1 && version.minor >= 1) || (version.major > 1)
        num_temp_sensor_channels = read(fid, Int16)
    end
    header.num_temp_sensor_channels = num_temp_sensor_channels
    
    # If data file is from GUI v1.3 or later, load eval board mode
    eval_board_mode = 0
    if (version.major == 1 && version.minor >= 3) || (version.major > 1)
        eval_board_mode = read(fid, Int16)
    end
    header.eval_board_mode = eval_board_mode
    
    # If data file is from v2.0 or later (Intan Recording Controller), load name of digital reference channel
    reference_channel = ""
    if version.major > 1
        reference_channel = readQString(fid)
    end
    header.reference_channel = reference_channel
    
    # If data file is from v2.0 or later (Intan Recording Controller), 128 samples in each data block. Otherwise, 60
    num_samples_per_data_block = 60
    if version.major > 1
        num_samples_per_data_block = 128
    end
    header.num_samples_per_data_block = num_samples_per_data_block
    
    # Place frequency-related information in data structure
    header.frequency_parameters.amplifier_sample_rate = header.sample_rate
    header.frequency_parameters.aux_input_sample_rate = header.sample_rate / 4
    header.frequency_parameters.supply_voltage_sample_rate = header.sample_rate / num_samples_per_data_block
    header.frequency_parameters.board_adc_sample_rate = header.sample_rate
    header.frequency_parameters.board_dig_in_sample_rate = header.sample_rate
    
    header.spike_triggers = SpikeTriggerStruct[]
    
    header.amplifier_channels = ChannelStruct[]
    header.aux_input_channels = ChannelStruct[]
    header.supply_voltage_channels = ChannelStruct[]
    header.board_adc_channels = ChannelStruct[]
    header.board_dig_in_channels = ChannelStruct[]
    header.board_dig_out_channels = ChannelStruct[]
    
    amplifier_index = 1
    aux_input_index = 1
    supply_voltage_index = 1
    board_adc_index = 1
    board_dig_in_index = 1
    board_dig_out_index = 1
    
    # Read signal summary from data file header
    number_of_signal_groups = read(fid, Int16)
    
    for signal_group = 1:number_of_signal_groups
        
        signal_group_name = readQString(fid)
        signal_group_prefix = readQString(fid)
        signal_group_enabled = read(fid, Int16)
        signal_group_num_channels = read(fid, Int16)
        signal_group_num_amp_channels = read(fid, Int16)
        
        if (signal_group_num_channels > 0) && (signal_group_enabled > 0)
            for signal_channel = 1:signal_group_num_channels
                new_trigger_channel = SpikeTriggerStruct(0, 0, 0, 0)
                
                new_channel = ChannelStruct("", "", 0, 0, 0, 0, "", "", 0, 0.0, 0.0)
                
                new_channel.port_name = signal_group_name
                new_channel.port_prefix = signal_group_prefix
                new_channel.port_number = signal_group
                
                new_channel.native_channel_name = readQString(fid)
                new_channel.custom_channel_name = readQString(fid)
                new_channel.native_order = read(fid, Int16)
                new_channel.custom_order = read(fid, Int16)
                signal_type = read(fid, Int16)
                channel_enabled = read(fid, Int16)
                new_channel.chip_channel = read(fid, Int16)
                new_channel.board_stream = read(fid, Int16)
                new_trigger_channel.voltage_trigger_mode = read(fid, Int16)
                new_trigger_channel.voltage_threshold = read(fid, Int16)
                new_trigger_channel.digital_trigger_channel = read(fid, Int16)
                new_trigger_channel.digital_edge_polarity = read(fid, Int16)
                new_channel.electrode_impedance_magnitude = read(fid, Float32)
                new_channel.electrode_impedance_phase = read(fid, Float32)
                
                if channel_enabled > 0
                    if signal_type == 0
                        push!(header.amplifier_channels, new_channel)
                        push!(header.spike_triggers, new_trigger_channel)
                        amplifier_index = amplifier_index + 1
                    elseif signal_type == 1
                        push!(header.aux_input_channels, new_channel)
                        aux_input_index = aux_input_index + 1
                    elseif signal_type == 2
                        push!(header.supply_voltage_channels, new_channel)
                        supply_voltage_index = supply_voltage_index + 1
                    elseif signal_type == 3
                        push!(header.board_adc_channels, new_channel)
                        board_adc_index = board_adc_index + 1
                    elseif signal_type == 4
                        push!(header.board_dig_in_channels, new_channel)
                        board_dig_in_index = board_dig_in_index + 1
                    elseif signal_type == 5
                        push!(header.board_dig_out_channels, new_channel)
                        board_dig_out_index = board_dig_out_index + 1
                    else
                        error("Unknown channel type")
                    end
                end
            end
        end
    end
    
    # Summarize contents of data file
    header.num_amplifier_channels = amplifier_index - 1
    header.num_aux_input_channels = aux_input_index - 1
    header.num_supply_voltage_channels = supply_voltage_index - 1
    header.num_board_adc_channels = board_adc_index - 1
    header.num_board_dig_in_channels = board_dig_in_index - 1
    header.num_board_dig_out_channels = board_dig_out_index - 1
    
    return header
end

function data_to_result(header, data, data_present)
    # Moves the header and data (if present) into a common object
    result = Result()
    result.t_amplifier = data.t_amplifier
    result.t_aux_input = data.t_aux_input
    result.t_supply_voltage = data.t_supply_voltage
    result.t_board_adc = data.t_board_adc
    result.t_dig = data.t_dig
    result.t_temp_sensor = data.t_temp_sensor
    result.spike_triggers = header.spike_triggers
    
    result.notes = header.notes
    result.frequency_parameters = header.frequency_parameters
    
    result.reference_channel = header.reference_channel
    
    result.amplifier_channels = header.amplifier_channels
    result.amplifier_data = data.amplifier_data
    
    result.aux_input_channels = header.aux_input_channels
    result.aux_input_data = data.aux_input_data
    
    result.supply_voltage_channels = header.supply_voltage_channels
    result.supply_voltage_data = data.supply_voltage_data
    
    result.board_adc_channels = header.board_adc_channels
    result.board_adc_data = data.board_adc_data
    
    result.board_dig_in_channels = header.board_dig_in_channels
    result.board_dig_in_data = data.board_dig_in_data
    
    result.board_dig_out_channels = header.board_dig_out_channels
    result.board_dig_out_data = data.board_dig_out_data
    
    return result
end

function read_one_data_block(data, header, indices, fid)
    # Reads one 60 or 128 sample data block from fid into data, at the location indicated by indices
    
    # In version 1.2, we moved from saving timestamps as unsigned integers to signed integers to
    # accommodate negative (adjusted) timestamps for pretrigger data
    if (header.version.major == 1 && header.version.minor >= 2) || (header.version.major > 1)
        data.t_amplifier[indices.amplifier:(indices.amplifier + header.num_samples_per_data_block - 1)] = reinterpret(Int32, read(fid, header.num_samples_per_data_block * 4))
    else
        data.t_amplifier[indices.amplifier:(indices.amplifier + header.num_samples_per_data_block - 1)] = reinterpret(UInt32, read(fid, header.num_samples_per_data_block * 4))
    end
    
    if header.num_amplifier_channels > 0
        data.amplifier_data[:, indices.amplifier:(indices.amplifier + header.num_samples_per_data_block - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, header.num_samples_per_data_block * header.num_amplifier_channels * 2)), Int64(header.num_samples_per_data_block), Int64(header.num_amplifier_channels)))
    end
    
    if header.num_aux_input_channels > 0
        data.aux_input_data[:, indices.aux_input:(indices.aux_input + Int(header.num_samples_per_data_block / 4) - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, Int(header.num_samples_per_data_block * header.num_aux_input_channels * 2 / 4))), Int64(header.num_samples_per_data_block / 4), Int64(header.num_aux_input_channels)))
    end
    
    if header.num_supply_voltage_channels > 0
        data.supply_voltage_data[:, indices.supply_voltage] = reinterpret(UInt16, read(fid, header.num_supply_voltage_channels * 2))
    end
    
    if header.num_temp_sensor_channels > 0
        data.temp_sensor_data[:, indices.supply_voltage] = reinterpret(UInt16, read(fid, header.num_temp_sensor_channels * 2))
    end
    
    if header.num_board_adc_channels > 0
        data.board_adc_data[:, indices.board_adc:(indices.board_adc + header.num_samples_per_data_block - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, header.num_samples_per_data_block * header.num_board_adc_channels * 2)), Int64(header.num_samples_per_data_block), Int64(header.num_board_adc_channels)))
    end
    
    if header.num_board_dig_in_channels > 0
        data.board_dig_in_raw[indices.board_dig_in:(indices.board_dig_in + header.num_samples_per_data_block - 1)] = reinterpret(UInt16, read(fid, header.num_samples_per_data_block * 2))
    end
    
    if header.num_board_dig_out_channels > 0
        data.board_dig_out_raw[indices.board_dig_out:(indices.board_dig_out + header.num_samples_per_data_block - 1)] = reinterpret(UInt16, read(fid, header.num_samples_per_data_block * 2))
    end
end

function notch_filter(input, f_sample, f_notch, bandwidth)
    t_step = 1 / Float64(f_sample)
    f_c = f_notch * t_step

    l = length(input)

    # Calculate IIR filter parameters
    d = exp(-2 * pi * (bandwidth / 2) * t_step)
    b = (1 + d * d) * cos(2 * pi * f_c)
    a0 = 1
    a1 = -b
    a2 = d * d
    a = (1 + d * d) / 2
    b0 = 1
    b1 = -2 * cos(2 * pi * f_c)
    b2 = 1

    output = Vector{Float64}(undef, length(input))
    output[1] = input[1]
    output[2] = input[2]

    #= (If filtering a continuous data stream, change output[1] and output[2] to the previous final two values of out.) =#
    #= Run filter =#
    for k = 3:l
        output[k] = (a*b2*input[k-2] + a*b1*input[k-1] + a*b0*input[k] - a2*output[k-2] - a1*output[k-1])/a0
    end

    return output
end


# Define find_channel_in_group function
function find_channel_in_group(channel_name, signal_group)
    #for i in range(1, length(signal_group))
    for (count, this_channel) in enumerate(signal_group)
        if this_channel.custom_channel_name == channel_name
            return true, count
        end
    end
    return false, 0
end


# Define find_channel_in_header function
function find_channel_in_header(channel_name, header)
    # Look through all present signal groups
    
    # 1. Look through amplifier_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.amplifier_channels)
    if channel_found
        return true, "amplifier_channels", channel_index
    end
    
    # 2. Look through aux_input_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.aux_input_channels)
    if channel_found
        return true, "aux_input_channels", channel_index
    end
    
    # 3. Look through supply_voltage_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.supply_voltage_channels)
    if channel_found
        return true, "supply_voltage_channels", channel_index
    end
    
    # 4. Look through board_adc_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.board_adc_channels)
    if channel_found
        return true, "board_adc_channels", channel_index
    end
    
    # 5. Look through board_dig_in_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.board_dig_in_channels)
    if channel_found
        return true, "board_dig_in_channels", channel_index
    end
    
    # 6. Look through board_dig_out_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.board_dig_out_channels)
    if channel_found
        return true, "board_dig_out_channels", channel_index
    end
    
    return false, "", 0
    
end


# Define plot_channel function
function plot_channel(channel_name, result)
    # Find channel that corresponds to this name
    (channel_found, signal_type, signal_index) = find_channel_in_header(channel_name, result)
        
    # Plot this channel
    if channel_found
        plotly()
        
        if signal_type == "amplifier_channels"
            y_label = "Voltage (microVolts)"
            t_vector = result.t_amplifier
            data_vector = result.amplifier_data[signal_index, :]
        elseif signal_type == "aux_input_channels"
            y_label = "Voltage (Volts)"
            t_vector = result.t_aux_input
            data_vector = result.aux_input_data[signal_index, :]
        elseif signal_type == "supply_voltage_channels"
            y_label = "Voltage (Volts)"
            t_vector = result.t_supply_voltage
            data_vector = result.supply_voltage_data[signal_index, :]
        elseif signal_type == "board_adc_channels"
            y_label = "Voltage (Volts)"
            t_vector = result.t_board_adc
            data_vector = result.board_adc_data[signal_index, :]
        elseif signal_type == "board_dig_in_channels"
            y_label = "Digital In Events (High or Low)"
            t_vector = result.t_dig
            data_vector = result.board_dig_in_data[signal_index, :]
        elseif signal_type == "board_dig_out_channels"
            y_label = "Digital Out Events (High or Low)"
            t_vector = result.t_dig
            data_vector = result.board_dig_out_data[signal_index, :]
        else
            error("Plotting not possible; signal type ", signal_type, " not found")
        end
        
        display(plot(t_vector[:], data_vector[:], title = channel_name, xlabel = "Time (s)", ylabel = y_label, legend = false))
        
    else
        error("Plotting not possible; channel ", channel_name, " not found")
    end
end

# Define load_file function
function load_file(filename)
    # Start timing
    start = time()
    
    # Open file
    fid = open(filename, "r")
    filesize = stat(filename).size
    
    # Read file header
    header = read_header(fid)
    
    # Output a summary of recorded data
    println("Found ", header.num_amplifier_channels, " amplifier channel", plural(header.num_amplifier_channels))
    println("Found ", header.num_aux_input_channels, " auxiliary input channel", plural(header.num_aux_input_channels))
    println("Found ", header.num_supply_voltage_channels, " supply voltage channel", plural(header.num_supply_voltage_channels))
    println("Found ", header.num_board_adc_channels, " board ADC channel", plural(header.num_board_adc_channels))
    println("Found ", header.num_board_dig_in_channels, " board digital input channel", plural(header.num_board_dig_in_channels))
    println("Found ", header.num_board_dig_out_channels, " board digital output channel", plural(header.num_board_dig_out_channels))
    println("Found ", header.num_temp_sensor_channels, " temperature sensor channel", plural(header.num_temp_sensor_channels), "\n")
    
    # Determine how many samples the data file contains
    bytes_per_block = get_bytes_per_data_block(header)
    
    # Calculate how many data blocks are present
    data_present = 0
    bytes_remaining = filesize - position(fid)
    if bytes_remaining > 0
        data_present = 1
    end
    
    if bytes_remaining % bytes_per_block != 0
        error("Something is wrong with file size: should have a whole number of data blocks")
    end
    
    num_data_blocks = Int(bytes_remaining / bytes_per_block)
    
    # Calculate how many samples of each signal type are present
    num_amplifier_samples = Int(header.num_samples_per_data_block * num_data_blocks)
    num_aux_input_samples = Int((header.num_samples_per_data_block / 4) * num_data_blocks)
    num_supply_voltage_samples = Int(1 * num_data_blocks)
    num_board_adc_samples = Int(header.num_samples_per_data_block * num_data_blocks)
    num_board_dig_in_samples = Int(header.num_samples_per_data_block * num_data_blocks)
    num_board_dig_out_samples = Int(header.num_samples_per_data_block * num_data_blocks)
    
    # Calculate how much time has been recorded
    record_time = num_amplifier_samples / header.sample_rate
    
    # Output a summary of contents of header file
    if data_present > 0
        @printf("File contains %0.3f seconds of data. Amplifiers were sampled at %0.2f kS/s.\n", record_time, header.sample_rate / 1000)
    else
        @printf("Header file contains no data. Amplifiers were sampled at %0.2f kS/s.\n", header.sample_rate / 1000)
    end
    
    if data_present > 0
        # Pre-allocate memory for data
        println("Allocating memory for data...\n")
        
        data = Data()
        if (header.version.major == 1 && header.version.minor >= 2) || (header.version.major > 1)
            data.t_amplifier = zeros(Int32, 1, num_amplifier_samples)
        else
            data.t_amplifier = zeros(UInt32, 1, num_amplifier_samples)
        end
        
        data.amplifier_data = zeros(UInt16, header.num_amplifier_channels, num_amplifier_samples)
        data.aux_input_data = zeros(Float64, header.num_aux_input_channels, num_aux_input_samples)
        data.supply_voltage_data = zeros(Float64, header.num_supply_voltage_channels, num_supply_voltage_samples)
        data.temp_sensor_data = zeros(Float64, header.num_temp_sensor_channels, num_supply_voltage_samples)
        data.board_adc_data = zeros(Float64, header.num_board_adc_channels, num_board_adc_samples)
        data.board_dig_in_data = zeros(Int16, header.num_board_dig_in_channels, num_board_dig_in_samples)
        data.board_dig_in_raw = Vector{UInt16}(undef, num_board_dig_in_samples)
        data.board_dig_out_data = zeros(Int16, header.num_board_dig_out_channels, num_board_dig_out_samples)
        data.board_dig_out_raw = Vector{UInt16}(undef, num_board_dig_out_samples)
        
        # Read sampled data from file
        println("Reading data from file...")
        
        # Initialize indices used in looping
        indices = Indices()
        
        print_increment = 10
        percent_done = print_increment
        for i = 1:num_data_blocks
            read_one_data_block(data, header, indices, fid)
            
            # Increment indices
            indices.amplifier = header.num_samples_per_data_block + indices.amplifier
            indices.aux_input = Int(header.num_samples_per_data_block / 4) + indices.aux_input
            indices.supply_voltage = 1 + indices.supply_voltage
            indices.board_adc = header.num_samples_per_data_block + indices.board_adc
            indices.board_dig_in = header.num_samples_per_data_block + indices.board_dig_in
            indices.board_dig_out = header.num_samples_per_data_block + indices.board_dig_out
            
            fraction_done = 100 * (1.0 * i / num_data_blocks)
            if fraction_done >= percent_done
                println(percent_done, "% done...")
                percent_done = percent_done + print_increment
            end
            
        end
        
        # Make sure we have read exactly the right amount of data
        bytes_remaining = filesize - position(fid)
        if bytes_remaining != 0
            error("Error: End of file not reached.")
        end
    else
        data = nothing
    end
    
    # Close data file
    close(fid)
    
    if data_present > 0
        println("Parsing data...\n")
        
        # Extract digital input channels to separate variables
        for i = 1 : header.num_board_dig_in_channels
            mask = 2^header.board_dig_in_channels[i].native_order
            data.board_dig_in_data[i,:] = (x -> (x > 0 ? 1 : 0)).(data.board_dig_in_raw .& mask)
        end
        
        # Extract digital output channels to separate variables
        for i = 1 : header.num_board_dig_out_channels
            mask = 2^header.board_dig_out_channels[i].native_order
            data.board_dig_out_data[i,:] = (x -> (x > 0 ? 1 : 0)).(data.board_dig_out_raw .& mask)
        end
        
        # Scale voltage levels appropriately
        data.amplifier_data = 0.195 .* (data.amplifier_data .- 32768) # units = microvolts
        data.aux_input_data = 37.4e-6 .* data.aux_input_data # units = volts
        data.supply_voltage_data = 74.8e-6 .* data.supply_voltage_data # units = volts
        
        if header.eval_board_mode == 1
            data.board_adc_data = 152.59e-6 .* (data.board_adc_data .- 32768) # units = volts
        elseif header.eval_board_mode == 13
            data.board_adc_data = 312.5e-6 .* (data.board_adc_data .- 32768) # units = volts
        else
            data.board_adc_data = 50.354e-6 .* data.board_adc_data # units = volts
        end
        
        data.temp_sensor_data = 0.01 .* data.temp_sensor_data # units = deg C
        
        # Check for gaps in timestamps
        num_gaps = sum(diff(data.t_amplifier, dims=2)[1, :] .!= 1)
        if num_gaps == 0
            println("No missing timestamps in data.")
        else
            println("Warning: ", num_gaps, " gaps in timestamp data found. Time scale will not be uniform!")
        end
        
        # Scale time steps (units = seconds)
        data.t_amplifier = data.t_amplifier / header.sample_rate
        max_num_samples = length(data.t_amplifier)
        data.t_aux_input = Array{Float64,2}(undef, 1, Int(max_num_samples/4))
        data.t_aux_input[1, :] = data.t_amplifier[1, 1:4:max_num_samples]
        data.t_supply_voltage = Array{Float64,2}(undef, 1, Int(max_num_samples / header.num_samples_per_data_block))
        data.t_supply_voltage[1, :] = data.t_amplifier[1:header.num_samples_per_data_block:max_num_samples]
        data.t_board_adc = data.t_amplifier
        data.t_dig = data.t_amplifier
        data.t_temp_sensor = data.t_supply_voltage
        
        # If the software notch filter was selected during the recording, apply the
        # same notch filter to amplifier data here
        if header.frequency_parameters.notch_filter_frequency > 0 && header.version.major < 3
            println("Applying notch filter...")
            
            print_increment = 10
            percent_done = print_increment
            for i = 1 : header.num_amplifier_channels
                data.amplifier_data[i, :] = notch_filter(data.amplifier_data[i, :], header.sample_rate, header.frequency_parameters.notch_filter_frequency, 10)
                fraction_done = 100 * (i / header.num_amplifier_channels)
                if fraction_done >= percent_done
                    println(percent_done, "% done...")
                    percent_done = percent_done + print_increment
                end
            end
            
        end
        
    end
    
    # Move variables to result struct
    result = data_to_result(header, data, data_present)
    
    elapsed = time() - start
    @printf("Done! Elapsed time: %0.1f seconds\n", elapsed)
    
    return (result, data_present)
end


# Define print_all_channel_names function
function print_all_channel_names(result)
    # Print all amplifier_channels
    print_names_in_group(result.amplifier_channels)
    
    # Print all aux_input_channels
    print_names_in_group(result.aux_input_channels)
    
    # Print all supply_voltage_channels
    print_names_in_group(result.supply_voltage_channels)
    
    # Print all board_adc_channels
    print_names_in_group(result.board_adc_channels)
    
    # Print all board_dig_in_channels
    print_names_in_group(result.board_dig_in_channels)
    
    # Print all board_dig_out_channels
    print_names_in_group(result.board_dig_out_channels)
end


# Define function print_names_in_group
function print_names_in_group(signal_group)
    for this_channel in signal_group
        println(this_channel.custom_channel_name)
    end
end