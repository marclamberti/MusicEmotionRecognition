#include <algorithm>
#include <array>
#include <aubio/aubio.h>
#include <cstdint>
#include <cmath>
#include <ctgmath>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <complex>
#include <valarray>

#include "xtract/libxtract.h"
#include "xtract/xtract_scalar.h"
#include "xtract/xtract_helper.h"
#include "WaveFile.h"

#include "../tools/LibXtract-master/src/fft.h"

#include "../tools/LibXtract-master/xtract/libxtract.h"
#include "../tools/LibXtract-master/src/xtract_macros_private.h"
#include "../tools/LibXtract-master/src/xtract_globals_private.h"

namespace {

	const unsigned int kBlockSize = 1024;
	const unsigned int kNumberOfJazzSong = 85;
	const unsigned int kNumberOfClassicalSong = 65;

	using Complex = std::complex<double>;
	using CArray = std::valarray<Complex>;
 
	const double PI = 3.141592653589793238460;

	enum Format {
		MIDI,
		WAV,
		NUMBER_OF_FORMAT,
	};

	std::vector<std::string> kFormatToString = {
		".mid",
		".wav",
	};

	enum Genders {
		CLASSICAL,
		JAZZ,
		NUMBER_OF_GENDERS
	};

	enum Modes {
		MINOR,
		MAJOR,
		NUMBER_OF_MODES
	};

	enum Keys {
		C,
		DB,
		D,
		EB,
		E,
		F,
		GB,
		G,
		AB,
		A,
		BB,
		B,
		NUMBER_OF_KEYS
	};

	std::map<std::string, enum Keys> kStringToKey = {
		{"C", Keys::C},
		{"DB", Keys::DB},
		{"D", Keys::D},
		{"EB", Keys::EB},
		{"E", Keys::E},
		{"F", Keys::F},
		{"GB", Keys::GB},
		{"G", Keys::G},
		{"AB", Keys::AB},
		{"A", Keys::A},
		{"BB", Keys::BB},
		{"B", Keys::B}
	};

	int kModeAndKeyIndex = 10;
	int kTempoIndex = 6;

	enum Features {
		BPM, // OK
		AVERAGE_ENERGY, // OK
		ENERGY_STANDARD_DEVIATION, // OK
		AVERAGE_FUNDAMENTAL_FREQUENCY,
		FUNDAMENTAL_FREQUENCY_STANDARD_DEVIATION,
		NUMBER_OF_FREQUENCIES_HIGHER_THAN_AVEREAGE_FUNDAMENTAL_FREQUENCY,
		AVERAGE_CENTROID, // MIDOK
		CENTROID_STANDARD_DEVIATION, // MIDOK
		MODE, // OK
		KEY, // OK
		GENDER, // OK
		NUMBER_OF_FEATURES
	};

	const int kValenceIndex = 11;
	const int kArousalIndex = 12;
	const int kLabelIndex = std::min(kArousalIndex, kValenceIndex);

	enum LabelTypes {
		AROUSAL,
		VALENCE,
		MULTILABEL,
		NUMBER_OF_LABEL_TYPES,
	};

}

std::string	GetFileExtension(std::string file) {
	std::string::size_type idx = file.rfind('.');
	if (idx != std::string::npos)
		return file.substr(idx + 1);;
	return std::string("");
}

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<std::string> Split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::string 	ExecCommand(std::string const &cmd) {
	std::string output;
	
	FILE *lsofFile_p = popen(cmd.c_str(), "r");
  	if (lsofFile_p)
  	{  
  		char buffer[1024];
  		memset(buffer, 0, 1024);
  		// For now we skip the first line.
  		char *line_p = fgets(buffer, sizeof(buffer), lsofFile_p);
  		line_p = fgets(buffer, sizeof(buffer), lsofFile_p);
  		pclose(lsofFile_p);
  		return line_p;
  	}

	return output;
}

std::string	IdToSongName(int song_id, enum Format f) {
	std::string song_number;
	std::string song_genre;

	if (song_id <= kNumberOfClassicalSong) {
		song_genre = "clas";
		if (song_id >= 10)
			song_number = std::to_string(song_id);
		else
			song_number = "0" + std::to_string(song_id);
	} else {
		song_genre = "jazz";
		if (song_id - kNumberOfClassicalSong >= 10)
			song_number = std::to_string(song_id - kNumberOfClassicalSong);
		else
			song_number = "0" + std::to_string(song_id - kNumberOfClassicalSong);
	}
	return song_genre + "_" + song_number + kFormatToString[f];
}

std::string GetInfoFromMIDI(std::string const &song_path) {
	char buff[512];
	memset(buff, 0, 512);
	getcwd(buff, 512);

	std::string cmd = "metamidi -l ";
	std::string cwd = buff;
	return ExecCommand(cmd + cwd + "/" + song_path);
}

enum Genders ExtractGender(int song_id) {
	return song_id <= kNumberOfClassicalSong ? Genders::CLASSICAL : Genders::JAZZ;
}

float ExtractBpm(std::vector<std::string> const &midi_info) {
	return std::atof(midi_info[kTempoIndex].c_str());
}

enum Modes ExtractMode(std::vector<std::string> const &midi_info) {
	if (midi_info[kModeAndKeyIndex].find('m') != std::string::npos)
		return Modes::MINOR;
	return Modes::MAJOR;
}

enum Keys ExtractKey(std::vector<std::string> const &midi_info) {
	std::string info_line = midi_info[kModeAndKeyIndex];
	std::transform(info_line.begin(), info_line.end(), info_line.begin(), ::toupper);
	std::vector<std::string> key_and_mode = Split(info_line, 'M');

	for (auto map_element : kStringToKey) {
		if (key_and_mode[0] == map_element.first)
			return map_element.second;
	}
	return Keys::NUMBER_OF_KEYS;
}

void	MidiFeatures(std::vector<float> &tuple, int song_id, std::string const &song_midi_path) {
	std::vector<std::string> info_in_midi = Split(GetInfoFromMIDI(song_midi_path), ';');
	tuple[BPM] = ExtractBpm(info_in_midi);
	tuple[MODE] = static_cast<float>(ExtractMode(info_in_midi));
	tuple[KEY] = static_cast<float>(ExtractKey(info_in_midi));
	tuple[GENDER] = static_cast<float>(ExtractGender(song_id));
}

void	ExtractFundamentalFrequencies(std::vector<float> &tuple, std::vector<double> &wav_data, std::uint64_t sample_rate) {
	std::vector<double>	fundamental_frequencies;
	double f0 = 0.;

}

// Cooley–Tukey FFT (in-place)
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];
 
    // conquer
    fft(even);
    fft(odd);
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

// Hamming window
std::vector<double> hamming(const std::vector<double> &data, int start, int len)
{
    std::vector<double>	window;
    for (int i = 0; i < len; ++i) {
    	double h = 0.54 - 0.46 * cos(2.0 * PI * i / (len - 1));
    	double result = h * data[start + i];
    	window.push_back(result);
    }
    return window;
}

// Spectrum centroid (Too slow)
// real	0m12.776s
// user	0m12.640s
// sys	0m0.106s
void	ExtractSpectrumCentroids(std::vector<float> &tuple, std::vector<double> &wav_data, std::uint64_t sample_rate) {
	std::vector<double> windowed;
	std::vector<double>	centroids;
	double sum_centroids = 0;

	std::cout << "Number of frames: " << wav_data.size() / kBlockSize << std::endl;
	for (int i = 0; (i + kBlockSize) < wav_data.size(); i += kBlockSize >> 1) {
	   	windowed = hamming(wav_data, i, kBlockSize);
	   	CArray x(windowed.size());
	   	for (int i = 0; i < windowed.size(); ++i) {
	   		x[i] = windowed[i];
	   	}
	   	fft(x);

	   	// Size between each bin
	   	//double bin_size = sample_rate / (double)kBlockSize;
	   	double centroid = 0.;
	   	double sum_frequency_and_magnitudes = 0.;
	   	double sum_magnitudes = 0.;
	   	// Get the frequency and the magnitude of each bin to compute the spectral centroid
	   	for (int bin = 0; bin < windowed.size() / 2; ++bin) {
	   		double frequency = (double)bin * (double)sample_rate / (double)windowed.size();
	   		double real = std::real(x[bin]);
	   		double imag = std::imag(x[bin]);
	   		double magnitude = std::sqrt(std::pow(real, 2) + std::pow(imag, 2));
	   		sum_frequency_and_magnitudes += frequency * magnitude;
	   		sum_magnitudes += magnitude;
	   		//std::cout << "bin: " << bin << " real: " << real << " imag: " << imag << " magnitude: " << magnitude << " phase (the angle between the two): " << phase << std::endl;
	   		//std::cout << "frequency: " << frequency << std::endl;
	   	}
	   	centroid = sum_frequency_and_magnitudes / sum_magnitudes;
	   	centroids.push_back(centroid);
	   	sum_centroids += centroid;
	   	//std::cout << "centroid : " << centroid << std::endl;
   }

 	// Compute the mean
 	double mean = sum_centroids / centroids.size();
 	tuple[AVERAGE_CENTROID] = static_cast<float>(mean);

 	// Compute standard deviation
 	double sum = 0.;
 	double std_deviation = 0.;
 	for (double centroid : centroids) {
 		sum += std::pow(centroid - mean, 2);
 	}
 	std_deviation = std::sqrt(sum / centroids.size());
 	tuple[CENTROID_STANDARD_DEVIATION] = static_cast<float>(std_deviation);
}

std::vector<double> FindEnergyInSamples(std::vector<double> &wav_file) {
	std::vector<double> energies;

	for (auto const &sample : wav_file) {
		double energy = pow(sample, 2);
		energies.push_back(energy);
	}
	return energies;
}

void	ExtractEnergy(std::vector<float> &tuple, std::vector<double> &wav_data) {
	std::vector<double> energy_vector = FindEnergyInSamples(wav_data);
	std::vector<double> means;

	auto energy_iterator = energy_vector.begin();
	for (uint64_t i = 0; i < wav_data.size(); i += kBlockSize) {
		unsigned int diff = wav_data.size() - i;
		int number_of_elements_available = std::min(diff, kBlockSize);
		double sum = std::accumulate(energy_iterator, energy_iterator + number_of_elements_available, 0.0);
		energy_iterator += kBlockSize;
		means.push_back(sum / number_of_elements_available);
	}

	double mean_of_means = std::accumulate(means.begin(), means.end(), 0.0) / means.size();
	double accum = 0.0;

	for (auto const &mean : means) {
		accum += (mean - mean_of_means) * (mean - mean_of_means);
	}
	double stdev = sqrt(accum / energy_vector.size());
	tuple[AVERAGE_ENERGY] = mean_of_means;
	tuple[ENERGY_STANDARD_DEVIATION] = static_cast<float>(stdev);
}

void	WavFeatures(std::vector<float> &tuple, std::string const &song_wav_path) {
	WaveFile wav_file(song_wav_path);
	if (!wav_file.IsLoaded()) {
		return;
	}

	float *data = (float *)wav_file.GetData();
	std::size_t bytes = wav_file.GetDataSize();
	std::uint64_t samples = bytes / sizeof(float);
	std::vector<double> wav_data(samples);
	std::uint64_t sample_rate = wav_file.GetSampleRate();
    for (std::uint64_t i = 0; i < samples; ++i) {
        wav_data[i] = (double)data[i];
    }
	//std::cout << "[WAVE] Bytes : " << bytes << std::endl;
	//std::cout << "[WAVE] Samples : " << samples << std::endl;
	//std::cout << "[WAVE] Sample Rate : " << sample_rate << std::endl;

	ExtractEnergy(tuple, wav_data);
	ExtractSpectrumCentroids(tuple, wav_data, sample_rate);
}

void	FillFeatures(std::vector<float> &tuple, int song_id, std::string const &song_midi_path, std::string const &song_wav_path) {
	//std::cout << "Name is : " << song_midi_path << std::endl;
	tuple.reserve(NUMBER_OF_FEATURES);
	//MidiFeatures(tuple, song_id, song_midi_path);
	WavFeatures(tuple, song_wav_path);
	//std::cout << "The beat is : " << tuple[BPM] << std::endl;
	//std::cout << "The mode is : " << tuple[MODE] << std::endl;
	//std::cout << "The key is : " << tuple[KEY] << std::endl;
	//std::cout << "The gender is : " << tuple[GENDER] << std::endl;
	//std::cout << "The average centroids is : " << tuple[AVERAGE_CENTROID] << std::endl;
	//std::cout << "The standard deviation of centroids is : " << tuple[CENTROID_STANDARD_DEVIATION] << std::endl;
	//std::cout << "The average energy is : " << tuple[AVERAGE_ENERGY] << std::endl;
	//std::cout << "The standard deviation of energy is : " << tuple[ENERGY_STANDARD_DEVIATION] << std::endl;
}

std::string FormatTuple(std::vector<float> const &tuple, enum LabelTypes lt) {
	std::string formatted_tuple;

	if (lt == LabelTypes::AROUSAL || lt == LabelTypes::MULTILABEL)
		formatted_tuple += std::to_string(tuple[kArousalIndex]);
	if (lt == LabelTypes::VALENCE || lt == LabelTypes::MULTILABEL)
		formatted_tuple += std::to_string(tuple[kArousalIndex]);

	int i = 1;
	for (auto const &attribute : tuple) {
		if (i <= kLabelIndex)
			formatted_tuple += " " + std::to_string(i++) + ':' + std::to_string(attribute);
	}
	return formatted_tuple;
}

void FormatDatasetAndWriteInFile(std::vector<std::vector<float>> const &data_set, std::string const &filename,
								 enum LabelTypes lt) {
	std::ofstream output_file_stream(filename.c_str(), std::fstream::out);

	if (!output_file_stream.is_open()) {
		std::cerr << "Cannot open the file " << filename << std::endl;
		return;
	}

	for (auto const &tuple : data_set) {
		std::string line = FormatTuple(tuple, lt);
		std::cout << "The line : " << line << std::endl;
		line += '\n';
		output_file_stream << line;
	}
	output_file_stream.close();
}

int main(int ac, char **av){
	std::vector<std::vector<float>> data_set(kNumberOfClassicalSong + kNumberOfJazzSong); 
	
	if (ac < 2) {
		std::cout << "Usage:" << av[0] << " <training_directory>" << std::endl;
		return 1;
	}
	int song_id = 1;
	for (auto &tuple : data_set) {
		std::string song_midi_path = av[1] + std::string("/") + IdToSongName(song_id, Format::MIDI);
		std::string song_wav_path = av[1] + std::string("/") + IdToSongName(song_id, Format::WAV);
		if (song_id < 2)
			FillFeatures(tuple, song_id++, song_midi_path, song_wav_path);
	}
	//FormatDatasetAndWriteInFile(data_set, "test.txt", LabelTypes::AROUSAL);
    return 0;
}