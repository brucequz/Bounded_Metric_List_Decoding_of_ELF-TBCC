#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class OutputFile {
public:
    explicit OutputFile(const std::string& filename) {
        m_file.open(filename, std::ios::out);
        if (!m_file.is_open()) {
            throw std::runtime_error("Could not open file: " + filename);
        }
    }

    ~OutputFile() {
        if (m_file.is_open()) m_file.close();
    }

    // Function to write 2D vectors (Matrix/Table)
    template <typename T>
    void write2DVector(const std::vector<std::vector<T>>& data, char delimiter = ',') {
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                m_file << row[i];
                
                // Add delimiter unless it's the last element in the row
                if (i < row.size() - 1) {
                    m_file << delimiter;
                }
            }
            m_file << "\n"; // New line after each row
        }
    }

    // Disable copying
    OutputFile(const OutputFile&) = delete;
    OutputFile& operator=(const OutputFile&) = delete;

private:
    std::ofstream m_file;
};