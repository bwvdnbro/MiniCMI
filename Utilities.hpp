/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file Utilities.hpp
 *
 * @brief General functions that are used throughout the program and are not
 * part of the standard library.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <sstream>
#include <vector>

/**
 * @brief Utility functions that are not really related to a single class.
 */
namespace Utilities {
/**
 * @brief Return a std::vector that contains the indices that would sort the
 * given vector.
 *
 * @param values std::vector to sort.
 * @return std::vector containing the indices of the elements in the given
 * std::vector in an order that would sort the std::vector.
 */
template <typename _datatype_>
std::vector<uint_fast32_t> argsort(const std::vector<_datatype_> &values) {

  std::vector<uint_fast32_t> idx(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    idx[i] = i;
  }
  std::sort(idx.begin(), idx.end(), [&values](size_t i1, size_t i2) {
    return values[i1] < values[i2];
  });
  return idx;
}

/**
 * @brief Compose a filename made up by the given prefix and counter value,
 * appropriately zero padded.
 *
 * @param folder Folder to add to the filename.
 * @param prefix Prefix for the filename.
 * @param extension Extension for the filename.
 * @param counter Value of the counter.
 * @param padding Number of digits the counter should have.
 * @return std::string with format: "<prefix>XX<counter>XX.<extension>", where
 * the number of Xs is equal to padding.
 */
inline std::string compose_filename(const std::string &folder,
                                    const std::string &prefix,
                                    const std::string &extension,
                                    uint_fast32_t counter,
                                    uint_fast32_t padding) {
  std::stringstream namestring;
  if (!folder.empty()) {
    namestring << folder << "/";
  }
  namestring << prefix;
  namestring.fill('0');
  namestring.width(padding);
  namestring << counter;
  namestring << "." << extension;
  return namestring.str();
}
} // namespace Utilities

#endif // UTILITIES_HPP
