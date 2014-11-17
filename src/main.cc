#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <ctgmath>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <set>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <valarray>
#include <vector>
#include <random>

#include "WaveFile.h"

namespace {

	const unsigned int kBlockSize = 1024;

	using Complex = std::complex<double>;
	using CArray = std::valarray<Complex>;
 
	enum Format {
		MIDI,
		WAV,
		NUMBER_OF_FORMAT,
	};

	std::vector<std::string> kFormatToString = {
		".mid",
		"_trimmed.wav",
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
		BPM, // 0
		AVERAGE_ENERGY, // 1
		ENERGY_STANDARD_DEVIATION, // 2
		AVERAGE_FUNDAMENTAL_FREQUENCY, // 3
		FUNDAMENTAL_FREQUENCY_STANDARD_DEVIATION, // 4
		NUMBER_OF_FREQUENCIES_HIGHER_THAN_AVEREAGE_FUNDAMENTAL_FREQUENCY, // 5
		AVERAGE_CENTROID, // 6
		CENTROID_STANDARD_DEVIATION, // 7
		MODE, // 8
		KEY, // 9
		GENDER, // 10
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

	std::string features_to_string[] = {
		"BPM",
		"AVERAGE_ENERGY",
		"ENERGY_STANDARD_DEVIATION",
		"AVERAGE_FUNDAMENTAL_FREQUENCY",
		"FUNDAMENTAL_FREQUENCY_STANDARD_DEVIATION",
		"NUMBER_OF_FREQUENCIES_HIGHER_THAN_AVEREAGE_FUNDAMENTAL_FREQUENCY",
		"AVERAGE_CENTROID",
		"CENTROID_STANDARD_DEVIATION",
		"MODE",
		"KEY",
		"GENDER"
	};

	std::string kFinalDataSetDirectory = "final_datasets";
	std::string kTrainingSetExtension = ".dataset";
	std::string kTestSetExtension = ".test";
	std::string kCombinaisonsDirectory = "combinaisons";

	const char *labels_to_string[] = {
		"AROUSAL",
		"VALENCE",
		"MULTILABEL"
	};

	enum SetType {
		TRAINING_SET,
		TEST_SET,
		NUMBER_OF_SET_TYPE
	};

	using Tuple = std::vector<float>;
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

std::vector<std::string> 	ExecCommand(std::string const &cmd) {
	std::vector<std::string> output;
	
	FILE *lsofFile_p = popen(cmd.c_str(), "r");
  	if (lsofFile_p)
  	{  
  		char buffer[1024];
  		memset(buffer, 0, 1024);
  		// For now we skip the first line.
  	  	char *line_p;
  		while ((line_p = fgets(buffer, sizeof(buffer), lsofFile_p))) {
  			output.push_back(line_p);
  		}

  		line_p = fgets(buffer, sizeof(buffer), lsofFile_p);
  		pclose(lsofFile_p);
  		return output;
  	}
	return output;
}

void GetSongsInDirectory(const std::string &directory_path, std::set<std::string> &files) {
	DIR *dir;
	struct dirent *file;
	struct stat	st;

	if ((dir = opendir(directory_path.c_str())) != NULL) {
		while ((file = readdir(dir)) != NULL) {
			lstat(file->d_name, &st);
			if (S_ISREG(st.st_mode)) {
				std::string fn = file->d_name;
				if (fn.substr(fn.find_last_of(".") + 1) == "wav" || fn.substr(fn.find_last_of(".") + 1) == "mid") {
					files.insert(fn);
				}
			}
		}
		closedir(dir);
	} 
}

std::vector<std::string> GetInfoFromMIDI(std::string const &song_path) {
	char buff[512];
	memset(buff, 0, 512);
	getcwd(buff, 512);

	std::string cmd = "metamidi -l ";
	std::string cwd = buff;
	return ExecCommand(cmd + cwd + "/" + song_path);
}

enum Genders ExtractGender(const std::string &song_midi_path) {
	return (song_midi_path.find("jazz") != std::string::npos) ? Genders::JAZZ : Genders::CLASSICAL;
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

void	MidiFeatures(std::vector<float> &tuple, std::string const &song_midi_path) {
	std::vector<std::string> info_in_midi = Split(GetInfoFromMIDI(song_midi_path)[1], ';');
	tuple[BPM] = ExtractBpm(info_in_midi);
	tuple[MODE] = static_cast<float>(ExtractMode(info_in_midi));
	tuple[KEY] = static_cast<float>(ExtractKey(info_in_midi));
	tuple[GENDER] = static_cast<float>(ExtractGender(song_midi_path));
}

// Cooley–Tukey FFT (in-place)
void FFT(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];
 
    // conquer
    FFT(even);
    FFT(odd);
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

// Hamming window
void Hamming(const std::vector<float> &data, std::vector<float> &window, int start, int len)
{
    for (int i = 0; i < len; ++i) {
    	float h = 0.54 - 0.46 * cos(2.0 * M_PI * i / (len - 1));
    	float result = h * data[start + i];
    	window.push_back(result);
    }
}

/**
 *  Calculates exact frequency of component in f_bin.
 *  Assumes FFT has already been applied.
 */
float GetFrequenceFromBin(const CArray &fft_data, int f_bin, std::uint64_t sample_rate) {
	float f_phase;
    if (std::abs(std::complex<float>(std::real(fft_data[f_bin]), std::imag(fft_data[f_bin]))) < 1)
        return 0;

    std::complex<float> val = std::complex<float>(std::real(fft_data[f_bin]), std::imag(fft_data[f_bin]));
    f_phase = arg(val);

    float freq_per_bin = (float)sample_rate / (float)kBlockSize;
    float cf = f_bin * freq_per_bin;
    float phase_change = f_phase;
    float expected = cf * (float)(kBlockSize >> 1) / (float)sample_rate;

    float phase_diff = phase_change / (2.0 * M_PI) - expected;
    phase_diff -= floor(phase_diff);

    if ((phase_diff -= floor(phase_diff)) > 0.5)
        phase_diff -= 1;

    phase_diff *= 2 * M_PI;

    double freq_diff = phase_diff * freq_per_bin * ((float)kBlockSize / (float)(kBlockSize >> 1)) / (2 * M_PI);
    double freq = cf + freq_diff;
    return freq;
}

/**
 * If the input signal is a musical note, then its spectrum should consist of a series of peaks, 
 * corresponding to fundamental frequency with harmonic components at integer multiples of the fundamental frequency. 
 * Hence when we compress the spectrum a number of times (downsampling), and compare it with the original spectrum, 
 * we can see that the strongest harmonic peaks line up. The first peak in the original spectrum coincides with the 
 * second peak in the spectrum compressed by a factor of two, which coincides with the third peak in the spectrum 
 * compressed by a factor of three. Hence, when the various spectrums are multiplied together, the result will 
 * form clear peak at the fundamental frequency.
 */
float HPS(const CArray &fft_data, std::uint64_t sample_rate) {
    float max = 0;
    int f_bin = 0;

    //calculate max HPS - only covering 1/6 of framesize
    //downsampling by factor of 3 * 1/2 of framesize
    for (int i = 0; i < kBlockSize / 6; ++i) {
        int i2 = 2 * i;
        int i3 = 3 * i;
        float hps =
            std::abs(std::complex<float>(std::real(fft_data[i]), std::imag(fft_data[i]))) +
            0.8 * std::abs(std::complex<float>(std::real(fft_data[i2]), std::imag(fft_data[i2]))) +
            0.6 * std::abs(std::complex<float>(std::real(fft_data[i3]), std::imag(fft_data[i3])));
        if (max < hps) {
            max = hps;
            f_bin = i;
        }
    }
    return GetFrequenceFromBin(fft_data, f_bin, sample_rate);
}

void	ExtractSpectrumCentroids(std::vector<float> &tuple,
	const std::vector<float> &centroids, float sum_centroids) {

 	// Compute the mean
 	float mean = sum_centroids / centroids.size();
 	tuple[AVERAGE_CENTROID] = mean;

 	// Compute the standard deviation
 	float sum = 0.;
 	float std_deviation = 0.;
 	for (float centroid : centroids) {
 		sum += std::pow(centroid - mean, 2);
 	}
 	std_deviation = std::sqrt(sum / centroids.size());
 	tuple[CENTROID_STANDARD_DEVIATION] = std_deviation;
}

void	ExtractFundamentalFrequencies(std::vector<float> &tuple,
		const std::vector<float> &fundamental_frequencies, float sum_fundamental_frequencies)
{
	// Compute the mean
	float mean = sum_fundamental_frequencies / fundamental_frequencies.size();
	tuple[AVERAGE_FUNDAMENTAL_FREQUENCY] = mean;

	// Compute the standard deviation and the number of f0 higher than the mean
	float sum = 0.;
	float std_deviation = 0.;
	float num_f0_higher_than_mean = 0.;
	for (float f0 : fundamental_frequencies) {
		sum += std::pow(f0 - mean, 2);
		if (f0 > mean) {
			++num_f0_higher_than_mean;
		}
	}
	std_deviation = std::sqrt(sum / fundamental_frequencies.size());
	tuple[FUNDAMENTAL_FREQUENCY_STANDARD_DEVIATION] = std_deviation;
	tuple[NUMBER_OF_FREQUENCIES_HIGHER_THAN_AVEREAGE_FUNDAMENTAL_FREQUENCY] = num_f0_higher_than_mean;
}

// FFT (Too slow)
// real	0m12.776s
// user	0m12.640s
// sys	0m0.106s
void	ExtractFeaturesUsingFFT(std::vector<float> &tuple, std::vector<float> &wav_data, std::uint64_t sample_rate) {
	std::vector<float> 	windowed;
	std::vector<float>	fundamental_frequencies;
	std::vector<float>	centroids;
	float sum_centroids = 0.;
	float sum_fundamental_frequencies = 0.;

	std::cout << "Number of frames: " << wav_data.size() / kBlockSize << std::endl;
	for (int start_block = 0; (start_block + kBlockSize) < wav_data.size(); start_block += kBlockSize >> 1) {
		windowed.clear();

	   	Hamming(wav_data, windowed, start_block, kBlockSize);
	   	CArray fft_data(windowed.size());
	   	for (int i = 0; i < windowed.size(); ++i) {
	   		fft_data[i] = windowed[i];
	   	}
	   	FFT(fft_data);

	   	float centroid = 0.;
	   	float sum_frequency_and_magnitudes = 0.;
	   	float sum_magnitudes = 0.;
	   	// Get the frequency and the magnitude of each bin to compute the spectral centroid
	   	for (int bin = 0; bin < windowed.size() / 2; ++bin) {
	   		float frequency = (float)bin * (float)sample_rate / (float)windowed.size();
	   		float real = std::real(fft_data[bin]);
	   		float imag = std::imag(fft_data[bin]);
	   		float magnitude = std::sqrt(std::pow(real, 2) + std::pow(imag, 2));
	   		sum_frequency_and_magnitudes += frequency * magnitude;
	   		sum_magnitudes += magnitude;
	   	}
	   	// Get the fundamental frequency of the window
	   	float freq = HPS(fft_data, sample_rate);
	   	if (freq) {
	   		fundamental_frequencies.push_back(freq);
	   		sum_fundamental_frequencies += freq;
	   	}
	   	centroid = sum_frequency_and_magnitudes / sum_magnitudes;
	   	centroids.push_back(centroid);
	   	sum_centroids += centroid;
   }
   ExtractSpectrumCentroids(tuple, centroids, sum_centroids);
   ExtractFundamentalFrequencies(tuple, fundamental_frequencies, sum_fundamental_frequencies);
}

std::vector<float> FindEnergyInSamples(std::vector<float> &wav_file) {
	std::vector<float> energies;
	for (auto const &sample : wav_file) {
		float energy = pow(sample, 2);
		energies.push_back(energy);
	}
	return energies;
}

void	ExtractEnergy(std::vector<float> &tuple, std::vector<float> &wav_data) {
	std::vector<float> energy_vector = FindEnergyInSamples(wav_data);
	std::vector<float> means;

	auto energy_iterator = energy_vector.begin();
	for (uint64_t i = 0; i < wav_data.size(); i += kBlockSize / 2) {
		unsigned int diff = wav_data.size() - i;
		int number_of_elements_available = std::min(diff, kBlockSize);
		float sum = std::accumulate(energy_iterator, energy_iterator + number_of_elements_available, 0.0);
		energy_iterator += kBlockSize / 2;
		means.push_back(sqrt(sum / number_of_elements_available));
	}

	//for (auto mean : means) {
	//	std::cout << "m : " << mean << std::endl;
	//}
	float mean_of_means = std::accumulate(means.begin(), means.end(), 0.0) / means.size();
	float accum = 0.0;

	for (auto const &mean : means) {
		accum += (mean - mean_of_means) * (mean - mean_of_means);
	}
	float stdev = sqrt(accum / energy_vector.size());
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
	std::vector<float> wav_data(samples);
	std::uint64_t sample_rate = wav_file.GetSampleRate();
    for (std::uint64_t i = 0; i < samples; ++i) {
        wav_data[i] = (float)data[i];
    }
	//std::cout << "[WAVE] Bytes : " << bytes << std::endl;
	//std::cout << "[WAVE] Samples : " << samples << std::endl;
	//std::cout << "[WAVE] Sample Rate : " << sample_rate << std::endl;

	ExtractEnergy(tuple, wav_data);
	ExtractSpectrumCentroids(tuple, wav_data, sample_rate);
	ExtractFeaturesUsingFFT(tuple, wav_data, sample_rate);
}

void	FillFeatures(std::vector<float> &tuple, std::string const &song_midi_path,
					std::string const &song_wav_path) {
	//std::cout << "Name is : " << song_midi_path << std::endl;
	tuple.resize(NUMBER_OF_FEATURES + MULTILABEL);
	MidiFeatures(tuple, song_midi_path);
	WavFeatures(tuple, song_wav_path);
	// std::cout << "The beat is : " << tuple[BPM] << std::endl;
	// std::cout << "The mode is : " << tuple[MODE] << std::endl;
	// std::cout << "The key is : " << tuple[KEY] << std::endl;
	// std::cout << "The gender is : " << tuple[GENDER] << std::endl;
	// std::cout << "The average centroids is : " << tuple[AVERAGE_CENTROID] << std::endl;
	// std::cout << "The standard deviation of centroids is : " << tuple[CENTROID_STANDARD_DEVIATION] << std::endl;
	// std::cout << "The average energy is : " << tuple[AVERAGE_ENERGY] << std::endl;
	// std::cout << "The standard deviation of energy is : " << tuple[ENERGY_STANDARD_DEVIATION] << std::endl;
	// std::cout << "The average f0 is : " << tuple[AVERAGE_FUNDAMENTAL_FREQUENCY] << std::endl;
	// std::cout << "The standard deviation of f0 is : " << tuple[FUNDAMENTAL_FREQUENCY_STANDARD_DEVIATION] << std::endl;
	// std::cout << "The number of f0 higher than the average of f0 is : " << tuple[NUMBER_OF_FREQUENCIES_HIGHER_THAN_AVEREAGE_FUNDAMENTAL_FREQUENCY] << std::endl;
}

std::string FormatTuple(std::vector<float> const &tuple, enum LabelTypes lt, bool is_multiclass_data_set) {
	std::string formatted_tuple;

	const int kValenceOffsetFromEnd = 2;
	const int kArousalOffsetFromEnd = 1;

	if (is_multiclass_data_set) {

		if (lt == LabelTypes::AROUSAL)
			formatted_tuple += std::to_string(tuple[tuple.size() - kArousalOffsetFromEnd]);
		else if (lt == LabelTypes::VALENCE)
			formatted_tuple += std::to_string(tuple[tuple.size() - kValenceOffsetFromEnd]);
	}
	else {
		formatted_tuple += std::to_string(tuple[tuple.size() - 1]);
	}

	int i = 1;
	int number_of_label = is_multiclass_data_set ? 2 : 1;
	for (auto const &attribute : tuple) {
		if (i <= tuple.size() - number_of_label)
			formatted_tuple += " " + std::to_string(i++) + ':' + std::to_string(attribute);
	}
	return formatted_tuple;
}

int FindLastLabelIndex(std::vector<std::string> const &output) {
	int index = 0;
	for (auto const &line : output) {
		if  (line.find(':') != std::string::npos)
			return index - 1;
		++index;
	}
	return index;
}

void FillFeaturesFromFile(std::vector<std::string> const &output, std::vector<float> &tuple, int index) {
	int i = 0;
	for (std::string const &line : output) {
		if (i > index) {
			std::vector<std::string> splited_line = Split(line, ':');
			tuple.push_back(std::atof(splited_line[1].c_str()));
		}
		++i;
	}
}

void FillLabelsFromFile(std::vector<std::string> const &output, std::vector<float> &tuple, int index) {
	int i = 0;
	for (std::string const &line : output) {
		if (i <= index) {
			tuple.push_back(std::atof(line.c_str()));
		}
		++i;
	}
}

std::vector<std::vector<float>> CreateDataSetFromFile(std::string const &filename) {
	std::ifstream input_file(filename.c_str());
	std::vector<std::vector<float>> new_data_set;

	if (!input_file.is_open()) {
		std::cerr << "Cannot open the file " << filename << std::endl;
		return new_data_set;
	}

	std::string line;
	while (std::getline(input_file,line))
    {
    	//std::cout << line << std::endl;
    	Tuple tuple;
    	std::vector<std::string> output = Split(line, ' ');
    	int index = FindLastLabelIndex(output);
    	FillFeaturesFromFile(output, tuple, index);\
    	FillLabelsFromFile(output, tuple, index);
    	new_data_set.push_back(tuple);
    }
    input_file.close();
    return new_data_set;
}

void FormatDatasetAndWriteInFile(std::vector<std::vector<float>> const &data_set, std::string const &filename,
								 enum LabelTypes lt, bool is_multiclass_data_set) {
	std::ofstream output_file_stream(filename.c_str(), std::fstream::out);

	if (!output_file_stream.is_open()) {
		std::cerr << "Cannot open the file " << filename << std::endl;
		return;
	}

	for (auto const &tuple : data_set) {
		std::string line = FormatTuple(tuple, lt, is_multiclass_data_set);
		//std::cout << "The line : " << line << std::endl;
		line += '\n';
		output_file_stream << line;
	}
	output_file_stream.close();
}

std::vector<std::vector<float>> ScaleDataSet(std::vector<std::vector<float>> const &data_set, enum LabelTypes lt,
		bool is_multiclass_data_set) {
	std::string filename = "tmp_files/dataset.dataset";
	std::string output_filename = "tmp_files/scaled_dataset.dataset";

	//std::cout << "Write the datasets " << std::endl;
	FormatDatasetAndWriteInFile(data_set, filename, lt, is_multiclass_data_set);
	std::string cmd = "svm-scale " + filename + " > " + output_filename;
	//std::cout << "Run svm scaled the datasets " << std::endl;
	ExecCommand(cmd);
	return CreateDataSetFromFile(output_filename);
}


typedef struct s_classification_result {
	std::vector<Tuple> training_set;
	std::vector<Tuple> test_set;	
	enum LabelTypes lt;
	std::vector<int> features_ids;
	float accuracy;
} t_classification_result;

int Factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : Factorial(n - 1) * n;
}

void InitializeFeaturesInUse(std::vector<int> &features_in_use) {
	int i = features_in_use.size() - 1;
	for (auto &feature_id : features_in_use) {
		feature_id = i;
		--i;
	}
}

void UpdateFeaturesInUse(std::vector<int> &features_in_use) {
	for (int i = 0; i < features_in_use.size(); ++i) {
		if (i == 0)
			features_in_use[i] = features_in_use[i] + 1;

			int j;
			for (j = i; j < features_in_use.size() && 
				features_in_use[j] == static_cast<int>(NUMBER_OF_FEATURES) - j; ++j) {
				features_in_use[j + 1] += 1;
				features_in_use[j] = 0;
			}
			if (features_in_use.size() > 1) {
				for (j = (j == features_in_use.size()) ? j - 2 : j - 1; j >= 0; --j) {
					while (features_in_use[j] <= features_in_use[j + 1]) {
						features_in_use[j] += 1;
					}
				}
			}
	}
}

std::string FindFeatureFileName(std::vector<int> const &features_ids, enum SetType st) {
	std::string filename = "F" + std::to_string(features_ids.size());
	for (auto const &feature : features_ids) {
		filename += "_";
		filename += std::to_string(feature);//features_to_string[feature];
	}
	if (st == TRAINING_SET)
		filename += ".dataset";
	else if (st == TEST_SET)
		filename += ".test";
	return filename;
}

std::vector<std::vector<float>> CreateDataSetWithSelectedFeatures(std::vector<std::vector<float>> const &data_set,
		std::vector<int> features_in_use) {

	std::vector<std::vector<float>> new_data_set;
	for (auto const &tuple : data_set) {
		Tuple tuple_with_new_features;
		for (auto const &feature : features_in_use) {
			tuple_with_new_features.push_back(tuple[feature]);
		}
		//std::cout << "Le tuple a une taille de : " << tuple.size() << std::endl;
		//for (auto value : tuple) {
		//	std::cout << " " << value;
		//}
		// //std::cout << std::endl;
		// if (lt == VALENCE)
		// 	tuple_with_new_features.push_back(tuple.size() - 1);
		// else
		tuple_with_new_features.push_back(tuple[tuple.size() - 1]);
		new_data_set.push_back(tuple_with_new_features);
	}
	return new_data_set;
}


// Modify s_classification_result and add test set and training set.
// Add Test set and training set here.
std::vector<struct s_classification_result *> GenerateCombinaisons(
	std::vector<std::vector<float>> const &training_set,
	std::vector<std::vector<float>> const &test_set, int number_of_features_used, enum LabelTypes lt) {
	std::vector<struct s_classification_result *> combinaisons;

	int const number_of_combinaison = Factorial(NUMBER_OF_FEATURES) /
		(Factorial(number_of_features_used) * Factorial(NUMBER_OF_FEATURES - number_of_features_used));
	std::vector<int> features_in_use(number_of_features_used);
	InitializeFeaturesInUse(features_in_use);

	for (int i = 0; i < number_of_combinaison; ++i) {
		t_classification_result *result = new t_classification_result;
		result->lt = lt;
		result->features_ids = features_in_use;
		result->training_set = CreateDataSetWithSelectedFeatures(training_set, features_in_use);
		result->test_set = CreateDataSetWithSelectedFeatures(test_set, features_in_use);


		std::string sub_directory = labels_to_string[lt];
		std::string training_filename = FindFeatureFileName(result->features_ids, TRAINING_SET);
		std::string training_file_path = kCombinaisonsDirectory + "/" + sub_directory + "/" + training_filename;
		std::string test_filename = FindFeatureFileName(result->features_ids, TEST_SET);
		std::string test_file_path = kCombinaisonsDirectory + "/" + sub_directory + "/" + test_filename;

		FormatDatasetAndWriteInFile(result->training_set, training_file_path, lt, false);
		FormatDatasetAndWriteInFile(result->test_set, test_file_path, lt, false);
		std::string kModelDirectory = "model_files";

		std::string kSVMOutputFileDirectory = "svm_output_files";

		std::string first_cmd = "svm-train -s 4 -t 2 " + training_file_path + " "
			+ kModelDirectory + "/" + sub_directory + "/" + training_filename + ".model";
		ExecCommand(first_cmd);
		std::string second_cmd = "svm-predict " + test_file_path + " "
			+ kModelDirectory + "/" + sub_directory + "/" + training_filename + ".model" + " "
			+ kSVMOutputFileDirectory + "/" + sub_directory + "/" + training_filename.substr(0, training_filename.find('.'));
		std::vector<std::string> output = ExecCommand(second_cmd);
		std::string accuracy_in_string = Split(output[0], ' ')[4];
		result->accuracy = std::atof(accuracy_in_string.c_str());
		combinaisons.push_back(result);
		if (i + 1 != number_of_combinaison)
			UpdateFeaturesInUse(features_in_use);
	}
	return combinaisons;
}

struct s_classification_result *GenerateAllCombinaisonsAndChooseTheBest(
	std::vector<std::vector<float>> const &training_set,
	std::vector<std::vector<float>> const &test_set, enum LabelTypes lt) {
	struct s_classification_result *best_result = NULL;

	for (int i = 1; i <= NUMBER_OF_FEATURES; ++i) {
		std::vector<struct s_classification_result*> results =
			GenerateCombinaisons(training_set, test_set, i, lt);
		for (auto &result : results) {
			if (!best_result)
				best_result = result;
			else if (best_result->accuracy > result->accuracy)
				best_result = result;
		}
		for (auto &result : results) {
			if (result != best_result)
				delete result;
		}
	}
	return best_result;
}

void FillTrainingSetAndTestSet(std::vector<std::vector<float>> const &data_set,
		std::vector<std::vector<float>> &training_set,
		std::vector<std::vector<float>> &test_set, int number_of_classical_songs, int number_of_jazz_songs) {

	std::vector<std::vector<float>>	tmp_data_set(data_set);

	//std::cout << "Dataset size : " << data_set.size() << std::endl;
	int i = 0;
	int max_test_set = number_of_classical_songs / 4;
	int max_training_set = number_of_classical_songs - max_test_set;
	while (i < max_training_set) {
		training_set.push_back(data_set[i]);
		++i;
	}
	while (i < max_training_set + max_test_set) {
		test_set.push_back(data_set[i]);
		++i;
	}

	max_test_set = number_of_jazz_songs / 4;
	max_training_set =  number_of_jazz_songs - max_test_set;
	while (i < number_of_classical_songs + max_training_set) {
		training_set.push_back(data_set[i]);
		++i;
	}
	while (i < number_of_classical_songs + max_training_set + max_test_set) {
		test_set.push_back(data_set[i]);
		++i;
	}

	// check if all sound has been treated
	int total_set = training_set.size() + test_set.size();
	if (total_set != data_set.size()) {
		std::cerr << "Some data set are not either in training set or test set, total set with training and test: " << total_set << " and number of datasets : " << data_set.size() << std::endl;
	}
}

void GenerateFinalTrainingAndTestSet(std::vector<std::vector<float>> const &training_set,
		std::vector<std::vector<float>> const &test_set, enum LabelTypes lt) {

	struct s_classification_result *best_result =
		GenerateAllCombinaisonsAndChooseTheBest(training_set, test_set, lt);


	std::cout << " ------------ THE BEST RESULT ------------ " << std::endl;
	std::cout << "Feature used : ";
	for (int i : best_result->features_ids) {
		std::cout << features_to_string[i] << " "  << i << " ";
	}
	std::cout << std::endl;
	std::cout << "The mean square error is : " << best_result->accuracy << std::endl;
	std::cout << " ------------------------------------------ " << std::endl;

	std::string training_set_filename = FindFeatureFileName(best_result->features_ids, TRAINING_SET);
	std::string test_set_filename = FindFeatureFileName(best_result->features_ids, TEST_SET);

	std::string training_set_final_path = kFinalDataSetDirectory + "/"
		+ std::string(labels_to_string[lt]) + kTrainingSetExtension;

	std::string test_set_final_path = kFinalDataSetDirectory +  "/"
		+ std::string(labels_to_string[lt]) + kTestSetExtension;

	std::string sub_directory = labels_to_string[lt];

	std::string cmd = "cp " + kCombinaisonsDirectory + "/" + sub_directory + "/" + training_set_filename + " "
		+ training_set_final_path;
	ExecCommand(cmd);
	cmd = "cp " + kCombinaisonsDirectory +  "/" + sub_directory + "/" + test_set_filename + " "
		+ test_set_final_path;
	ExecCommand(cmd);


	std::string first_cmd = "svm-train -s 4 -t 2 " + training_set_final_path + " "
		+ training_set_final_path + ".model";
	ExecCommand(first_cmd);
	std::string second_cmd = "svm-predict " + test_set_final_path + " "
		+ training_set_final_path + ".model" + " "
		+ training_set_final_path + ".output";
	ExecCommand(second_cmd);
}

void FillLabels(std::vector<std::vector<float>> &data_set) {
	std::string cmd = "cat merge_results_script/Merged_Results.txt";
	std::vector<std::string> output = ExecCommand(cmd);
	//std::cout << "La size de l'ouput est de : " << output.size() << std::endl;
	//std::cout << "La size du data_set est de : " << data_set.size() << std::endl;
	int i = 0;
	//std::cout << "La size du data set est de : " << data_set.size() << std::endl;
	for (auto &tuple : data_set) {
		std::vector<std::string> splited_line = Split(output[i], '\t');
		tuple[kValenceIndex] = std::atof(splited_line[0].c_str());
		tuple[kArousalIndex] = std::atof(splited_line[1].c_str());

		//std::cout << "La valence du tuple est de  " << tuple[kValenceIndex] << std::endl;
		//std::cout << "La arousal du tuple est de  " << tuple[kArousalIndex] << std::endl;
		++i;
	}
}

void	UpdateNumberSongsByGender(const std::vector<float> &tuple, int &number_of_classical_songs, int &number_of_jazz_songs) {
	if (tuple[GENDER] == Genders::JAZZ) {
		++number_of_jazz_songs;
	} else {
		++number_of_classical_songs;
	}
}

int main(int ac, char **av) {
	std::set<std::string> files;
	std::vector<std::vector<float>> data_set;
	std::vector<std::vector<float>> arousal_training_set;
	std::vector<std::vector<float>> arousal_test_set;
	std::vector<std::vector<float>> valence_training_set;
	std::vector<std::vector<float>> valence_test_set;
	int	number_of_classical_songs = 0;
	int number_of_jazz_songs = 0;

	if (ac < 2) {
		std::cout << "Usage:" << av[0] << " <training_directory>" << std::endl;
		return 1;
	}
	
	//int song_limit = 20;
	std::string	directory_path = av[1] + std::string("/");
	GetSongsInDirectory(directory_path, files);
	std::cout << "FILES: " << files.size() << std::endl;

	int i = 0;
	for (auto it = files.begin(); it != files.end(); //&& i < song_limit;
		 ++it) {//, ++i) {
		std::vector<float> tuple;
		std::string midi_file = *it;
		std::string wav_file = *(++it);
		std::string song_midi_path = directory_path + midi_file;
		std::string song_wav_path = directory_path + wav_file;
		std::cout << "Sound " << song_midi_path << " , " << song_wav_path << " in progress..." << std::endl;
		FillFeatures(tuple, song_midi_path, song_wav_path);
		UpdateNumberSongsByGender(tuple, number_of_classical_songs, number_of_jazz_songs);
		data_set.push_back(tuple);
	}
	std::cout << "Number of classical songs: " << number_of_classical_songs << ", number of jazz songs: " << number_of_jazz_songs << std::endl;
	std::cout << "Feature filled now filling labels." << std::endl;
	FillLabels(data_set);

	std::cout << "Scaling the datasets " << std::endl;

	std::vector<std::vector<float>> arousal_data_set = ScaleDataSet(data_set, AROUSAL, true);
	std::vector<std::vector<float>> valence_data_set = ScaleDataSet(data_set, VALENCE, true);

	std::cout << "Scaling the datasets endl" << std::endl;

	// std::cout << "labels filled now creating training and test set." << std::endl;
	FillTrainingSetAndTestSet(arousal_data_set, arousal_training_set, arousal_test_set, number_of_classical_songs, number_of_jazz_songs);
	FillTrainingSetAndTestSet(valence_data_set, valence_training_set, valence_test_set, number_of_classical_songs, number_of_jazz_songs);

	std::cout << "Size of arousal training set: " << arousal_training_set.size() << ", size of arousal_test set: " << arousal_test_set.size() << std::endl;		
	std::cout << "Size of valence training set: " << valence_training_set.size() << ", size of valence_test set: " << valence_test_set.size() << std::endl;		

	// std::cout << "Training and test set done now finding best combinaisons." << std::endl;
	GenerateFinalTrainingAndTestSet(valence_training_set, valence_test_set, VALENCE);
	GenerateFinalTrainingAndTestSet(arousal_training_set, arousal_test_set, AROUSAL);
    return 0;
}