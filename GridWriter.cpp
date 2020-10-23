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
 * @file GridWriter.cpp
 *
 * @brief GridWriter implementation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#include "GridWriter.hpp"
#include "HDF5Tools.hpp"
#include "Utilities.hpp"

#include <vector>

/**
 * @brief Constructor.
 *
 * @param prefix Prefix for the name of the file to write.
 * @param output_folder Name of the folder where output files should be placed.
 * @param padding Number of digits used for the counter in the filenames.
 * @param compression Compress the HDF5 output?
 */
GridWriter::GridWriter(std::string prefix, std::string output_folder,
                       uint_fast8_t padding, const bool compression)
    : _output_folder(output_folder), _prefix(prefix), _padding(padding),
      _compression(compression) {

  // turn off default HDF5 error handling: we catch errors ourselves
  HDF5Tools::initialize();
}

/**
 * @brief Write a snapshot for a split grid with hydro.
 *
 * @param grid_creator Grid.
 * @param counter Counter value to add to the snapshot file name.
 * @param params ParameterFile containing the run parameters that should be
 * written to the file.
 * @param time Simulation time (in s).
 */
void GridWriter::write(DensitySubGridCreator<HydroDensitySubGrid> &grid_creator,
                       const uint_fast32_t counter, double time) {

  std::string filename = Utilities::compose_filename(_output_folder, _prefix,
                                                     "hdf5", counter, _padding);

  const Box<> box = grid_creator.get_box();

  HDF5Tools::HDF5File file =
      HDF5Tools::open_file(filename, HDF5Tools::HDF5FILEMODE_WRITE);

  // write header
  HDF5Tools::HDF5Group group = HDF5Tools::create_group(file, "Header");
  CoordinateVector<> boxsize = box.get_sides();
  HDF5Tools::write_attribute<CoordinateVector<>>(group, "BoxSize", boxsize);
  int32_t dimension = 3;
  HDF5Tools::write_attribute<int32_t>(group, "Dimension", dimension);
  std::vector<uint32_t> flag_entropy(6, 0);
  HDF5Tools::write_attribute<std::vector<uint32_t>>(group, "Flag_Entropy_ICs",
                                                    flag_entropy);
  std::vector<double> masstable(6, 0.);
  HDF5Tools::write_attribute<std::vector<double>>(group, "MassTable",
                                                  masstable);
  int32_t numfiles = 1;
  HDF5Tools::write_attribute<int32_t>(group, "NumFilesPerSnapshot", numfiles);
  const uint64_t number_of_cells = grid_creator.number_of_cells();
  std::vector<uint32_t> numpart(6, 0);
  numpart[0] = static_cast<uint32_t>(number_of_cells);
  std::vector<uint32_t> numpart_high(6, 0);
  numpart_high[0] = static_cast<uint32_t>(number_of_cells >> 32);
  HDF5Tools::write_attribute<std::vector<uint32_t>>(group, "NumPart_ThisFile",
                                                    numpart);
  HDF5Tools::write_attribute<std::vector<uint32_t>>(group, "NumPart_Total",
                                                    numpart);
  HDF5Tools::write_attribute<std::vector<uint32_t>>(
      group, "NumPart_Total_HighWord", numpart_high);
  HDF5Tools::write_attribute<double>(group, "Time", time);
  HDF5Tools::close_group(group);

  // THE CODE BELOW IS A HACK
  // In the actual CMacIonize code, the Parameters block contains a copy of the
  // parameters in the actual parameter file
  // We only include some of them here because they are read by analysis scripts
  // that we want to use.
  group = HDF5Tools::create_group(file, "Parameters");
  const CoordinateVector<int_fast32_t> num_subgrids =
      grid_creator.get_subgrid_layout();
  const CoordinateVector<int_fast32_t> num_cells_per_subgrid =
      grid_creator.get_subgrid_cell_layout();
  const CoordinateVector<int_fast32_t> num_cells(
      num_subgrids.x() * num_cells_per_subgrid.x(),
      num_subgrids.y() * num_cells_per_subgrid.y(),
      num_subgrids.z() * num_cells_per_subgrid.z());
  std::stringstream arraystr;
  arraystr << "[" << num_cells.x() << "," << num_cells.y() << ","
           << num_cells.z() << "]";
  std::string array = arraystr.str();
  HDF5Tools::write_attribute<std::string>(group, "DensityGrid:number of cells",
                                          array);
  arraystr.str("");
  arraystr << "[" << num_subgrids.x() << "," << num_subgrids.y() << ","
           << num_subgrids.z() << "]";
  array = arraystr.str();
  HDF5Tools::write_attribute<std::string>(
      group, "DensitySubGridCreator:number of subgrids", array);
  arraystr.str("");
  arraystr << "[" << box.get_anchor().x() << " m, " << box.get_anchor().y()
           << " m, " << box.get_anchor().z() << " m]";
  array = arraystr.str();
  HDF5Tools::write_attribute<std::string>(group, "SimulationBox:anchor", array);
  arraystr.str("");
  arraystr << "[" << box.get_sides().x() << " m, " << box.get_sides().y()
           << " m, " << box.get_sides().z() << " m]";
  array = arraystr.str();
  HDF5Tools::write_attribute<std::string>(group, "SimulationBox:sides", array);
  HDF5Tools::close_group(group);
  // END OF HACK

  // write runtime parameters
  group = HDF5Tools::create_group(file, "RuntimePars");
  // an uint_fast32_t does not necessarily have the expected 32-bit size, while
  // we really need a 32-bit variable to write to the file
  uint32_t uint32_iteration = counter;
  HDF5Tools::write_attribute<uint32_t>(group, "Iteration", uint32_iteration);
  HDF5Tools::close_group(group);

  // write units, we use SI units everywhere
  group = HDF5Tools::create_group(file, "Units");
  double unit_current_in_cgs = 1.;
  double unit_length_in_cgs = 100.;
  double unit_mass_in_cgs = 1000.;
  double unit_temperature_in_cgs = 1.;
  double unit_time_in_cgs = 1.;
  HDF5Tools::write_attribute<double>(group, "Unit current in cgs (U_I)",
                                     unit_current_in_cgs);
  HDF5Tools::write_attribute<double>(group, "Unit length in cgs (U_L)",
                                     unit_length_in_cgs);
  HDF5Tools::write_attribute<double>(group, "Unit mass in cgs (U_M)",
                                     unit_mass_in_cgs);
  HDF5Tools::write_attribute<double>(group, "Unit temperature in cgs (U_T)",
                                     unit_temperature_in_cgs);
  HDF5Tools::write_attribute<double>(group, "Unit time in cgs (U_t)",
                                     unit_time_in_cgs);
  HDF5Tools::close_group(group);

  // write particles
  // to limit memory usage, we first create all datasets, and then add the data
  // in small blocks
  group = HDF5Tools::create_group(file, "PartType0");
  // the code below differs from actual CMacIonize, where the desired output
  // data can be selected in the parameter file
  // in MiniCMI, we simply hard-code the desired output
  uint_fast32_t number_of_vector_props = 0;
  uint_fast32_t number_of_scalar_props = 0;
  HDF5Tools::create_dataset<CoordinateVector<>>(group, "Coordinates",
                                                numpart[0], _compression);
  ++number_of_vector_props;
  HDF5Tools::create_dataset<double>(group, "NumberDensity", numpart[0],
                                    _compression);
  ++number_of_scalar_props;
  /// TEST VARIABLES
  HDF5Tools::create_dataset<double>(group, "TestDensity", numpart[0],
                                    _compression);
  ++number_of_scalar_props;
  HDF5Tools::create_dataset<double>(group, "TestNeighbourDensitySum",
                                    numpart[0], _compression);
  ++number_of_scalar_props;

  const uint_fast32_t blocksize = 10000;
  uint_fast32_t block_offset = 0;
  for (auto gridit = grid_creator.begin();
       gridit != grid_creator.original_end(); ++gridit) {

    const uint_fast32_t numblock =
        (*gridit).get_number_of_cells() / blocksize +
        ((*gridit).get_number_of_cells() % blocksize > 0);
    for (uint_fast32_t iblock = 0; iblock < numblock; ++iblock) {
      const uint_fast32_t offset = iblock * blocksize;
      const uint_fast32_t upper_limit = std::min(
          offset + blocksize, uint_fast32_t((*gridit).get_number_of_cells()));
      const uint_fast32_t thisblocksize = upper_limit - offset;

      // the code below differs from CMacIonize, where the number of vector and
      // scalar properties is configurable
      // in MiniCMI, these numbers are hard-coded and should match the number of
      // hard-coded output datasets of each type
      std::vector<std::vector<CoordinateVector<>>> vector_props(
          number_of_vector_props,
          std::vector<CoordinateVector<>>(thisblocksize));
      std::vector<std::vector<double>> scalar_props(
          number_of_scalar_props, std::vector<double>(thisblocksize));

      size_t index = 0;
      for (auto cellit = (*gridit).hydro_begin() + offset;
           cellit != (*gridit).hydro_begin() + upper_limit; ++cellit) {
        uint_fast8_t vector_index = 0;
        uint_fast8_t scalar_index = 0;
        // FILL ARRAYS WITH VARIABLES (code differs from CMacIonize)
        vector_props[vector_index][index] =
            cellit.get_cell_midpoint() - box.get_anchor();
        ++vector_index;
        scalar_props[scalar_index][index] =
            cellit.get_ionization_variables().get_number_density();
        ++scalar_index;
        scalar_props[scalar_index][index] =
            cellit.get_hydro_variables().get_test_density();
        ++scalar_index;
        scalar_props[scalar_index][index] =
            cellit.get_hydro_variables().get_test_neighbour_density_sum();
        ++scalar_index;
        ++index;
      }

      uint_fast8_t vector_index = 0;
      uint_fast8_t scalar_index = 0;
      // WRITE ARRAYS TO FILE (code differs from CMacIonize)
      HDF5Tools::append_dataset<CoordinateVector<>>(group, "Coordinates",
                                                    block_offset + offset,
                                                    vector_props[vector_index]);
      ++vector_index;
      HDF5Tools::append_dataset<double>(group, "NumberDensity",
                                        block_offset + offset,
                                        scalar_props[scalar_index]);
      ++scalar_index;
      HDF5Tools::append_dataset<double>(group, "TestDensity",
                                        block_offset + offset,
                                        scalar_props[scalar_index]);
      ++scalar_index;
      HDF5Tools::append_dataset<double>(group, "TestNeighbourDensitySum",
                                        block_offset + offset,
                                        scalar_props[scalar_index]);
      ++scalar_index;
    }
    block_offset += (*gridit).get_number_of_cells();
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}
