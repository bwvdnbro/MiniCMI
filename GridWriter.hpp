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
 * @file GridWriter.hpp
 *
 * @brief HDF5-file writer for the distributed grid.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef GRIDWRITER_HPP
#define GRIDWRITER_HPP

#include "DensitySubGridCreator.hpp"
#include "HydroDensitySubGrid.hpp"

#include <string>

/**
 * @brief HDF5-file writer for the distributed grid.
 */
class GridWriter {
private:
  /*! @brief Folder where the snapshots are placed. */
  const std::string _output_folder;

  /*! @brief Prefix of the name for the file to write. */
  const std::string _prefix;

  /*! @brief Number of digits used for the counter in the filenames. */
  const uint_fast8_t _padding;

  /*! @brief Compress the HDF5 output? */
  const bool _compression;

public:
  GridWriter(
      std::string prefix, std::string output_folder = std::string("."),
      uint_fast8_t padding = 3,
      const bool compression = false);

  virtual void write(DensitySubGridCreator< HydroDensitySubGrid > &grid_creator,
                     const uint_fast32_t counter,
                     double time = 0.);
};

#endif // GRIDWRITER_HPP
