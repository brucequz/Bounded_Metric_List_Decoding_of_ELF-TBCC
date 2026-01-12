#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

class FileHandler {
public:
    explicit FileHandler(std::string filename) : m_filename(std::move(filename)) {}

    // --- Writing Function ---
    template <typename T>
    void write2DVector(const std::vector<std::vector<T>>& data, char delimiter = ',') {
        std::ofstream outFile(m_filename, std::ios::out);
        if (!outFile.is_open()) {
            throw std::runtime_error("Could not open file for writing: " + m_filename);
        }

        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                outFile << row[i];
                if (i < row.size() - 1) outFile << delimiter;
            }
            outFile << "\n";
        }
    }

    // --- New Reading Function ---
    template <typename T>
    std::vector<std::vector<T>> read2DVector(char delimiter = ',') {
        std::vector<std::vector<T>> data;
        std::ifstream inFile(m_filename);

        if (!inFile.is_open()) {
            throw std::runtime_error("Could not open file for reading: " + m_filename);
        }

        std::string line;
        while (std::getline(inFile, line)) {
            std::vector<T> row;
            std::stringstream ss(line);
            std::string cell;

            // Extract each element based on the delimiter
            while (std::getline(ss, cell, delimiter)) {
                // Convert string to the template type T using a stringstream
                std::stringstream convert(cell);
                T value;
                if (convert >> value) {
                    row.push_back(value);
                }
            }

            if (!row.empty()) {
                data.push_back(row);
            }
        }
        return data;
    }

private:
    std::string m_filename;
};