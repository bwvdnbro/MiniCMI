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

  // write parameters
  // NEEDS IMPLEMENTATION (we at least need the parameters for the subgrid
  // layout to be able to parse the output files correctly...)
  //  group = HDF5Tools::create_group(file, "Parameters");
  //  for (auto it = params.begin(); it != params.end(); ++it) {
  //    std::string key = it.get_key();
  //    std::string value = it.get_value();
  //    HDF5Tools::write_attribute< std::string >(group, key, value);
  //  }
  //  HDF5Tools::close_group(group);

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
  //  for (int_fast32_t property = 0; property < DENSITYGRIDFIELD_NUMBER;
  //       ++property) {
  //    if (fields.field_present(property)) {
  //      const std::string name = DensityGridWriterFields::get_name(property);
  //      if (DensityGridWriterFields::get_type(property) ==
  //          DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
  //        HDF5Tools::create_dataset< CoordinateVector<> >(group, name,
  //        numpart[0],
  //                                                        _compression);
  //      } else {
  //        if (DensityGridWriterFields::is_ion_property(property)) {
  //          for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES; ++ion) {
  //            if (fields.ion_present(property, ion)) {
  //              const std::string prop_name = name + get_ion_name(ion);
  //              HDF5Tools::create_dataset< double >(group, prop_name,
  //              numpart[0],
  //                                                  _compression);
  //            }
  //          }
  //        } else if (DensityGridWriterFields::is_heating_property(property)) {
  //          for (int_fast32_t heating = 0; heating < NUMBER_OF_HEATINGTERMS;
  //               ++heating) {
  //            if (fields.heatingterm_present(property, heating)) {
  //              const std::string prop_name = name + get_ion_name(heating);
  //              HDF5Tools::create_dataset< double >(group, prop_name,
  //              numpart[0],
  //                                                  _compression);
  //            }
  //          }
  //        } else {
  //          HDF5Tools::create_dataset< double >(group, name, numpart[0],
  //                                              _compression);
  //        }
  //      }
  //    }
  //  }
  // CREATE DATASETS

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

      //      std::vector< std::vector< CoordinateVector<> > > vector_props(
      //          fields.get_field_count(DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE),
      //          std::vector< CoordinateVector<> >(thisblocksize));
      //      std::vector< std::vector< double > > scalar_props(
      //          fields.get_field_count(DENSITYGRIDFIELDTYPE_SCALAR_DOUBLE),
      //          std::vector< double >(thisblocksize));
      // PUT IN THE CORRECT NUMBERS HERE
      std::vector<std::vector<CoordinateVector<>>> vector_props(
          0, std::vector<CoordinateVector<>>(thisblocksize));
      std::vector<std::vector<double>> scalar_props(
          0, std::vector<double>(thisblocksize));

      size_t index = 0;
      for (auto cellit = (*gridit).hydro_begin() + offset;
           cellit != (*gridit).hydro_begin() + upper_limit; ++cellit) {
        uint_fast8_t vector_index = 0;
        uint_fast8_t scalar_index = 0;
        //        for (int_fast32_t property = 0; property <
        //        DENSITYGRIDFIELD_NUMBER;
        //             ++property) {
        //          if (fields.field_present(property)) {
        //            if (DensityGridWriterFields::get_type(property) ==
        //                DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
        //              vector_props[vector_index][index] =
        //                  DensityGridWriterFields::get_vector_double_value(
        //                      property, cellit, box.get_anchor());
        //              ++vector_index;
        //            } else {
        //              if (DensityGridWriterFields::is_ion_property(property))
        //              {
        //                for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES;
        //                ++ion) {
        //                  if (fields.ion_present(property, ion)) {
        //                    scalar_props[scalar_index][index] =
        //                        DensityGridWriterFields::get_scalar_double_ion_value(
        //                            property, ion, cellit);
        //                    ++scalar_index;
        //                  }
        //                }
        //              } else if (DensityGridWriterFields::is_heating_property(
        //                             property)) {
        //                for (int_fast32_t heating = 0; heating <
        //                NUMBER_OF_HEATINGTERMS;
        //                     ++heating) {
        //                  if (fields.heatingterm_present(property, heating)) {
        //                    scalar_props[scalar_index][index] =
        //                        DensityGridWriterFields::
        //                            get_scalar_double_heating_value(property,
        //                            heating,
        //                                                            cellit);
        //                    ++scalar_index;
        //                  }
        //                }
        //              } else {
        //                scalar_props[scalar_index][index] =
        //                    DensityGridWriterFields::get_scalar_double_value(property,
        //                                                                     cellit);
        //                ++scalar_index;
        //              }
        //            }
        //          }
        //        }
        // FILL ARRAYS WITH VARIABLES
        ++index;
      }

      uint_fast8_t vector_index = 0;
      uint_fast8_t scalar_index = 0;
      //      for (int_fast32_t property = 0; property <
      //      DENSITYGRIDFIELD_NUMBER;
      //           ++property) {
      //        if (fields.field_present(property)) {
      //          const std::string name =
      //          DensityGridWriterFields::get_name(property); if
      //          (DensityGridWriterFields::get_type(property) ==
      //              DENSITYGRIDFIELDTYPE_VECTOR_DOUBLE) {
      //            HDF5Tools::append_dataset< CoordinateVector<> >(
      //                group, name, block_offset + offset,
      //                vector_props[vector_index]);
      //            ++vector_index;
      //          } else {
      //            if (DensityGridWriterFields::is_ion_property(property)) {
      //              for (int_fast32_t ion = 0; ion < NUMBER_OF_IONNAMES;
      //              ++ion) {
      //                if (fields.ion_present(property, ion)) {
      //                  const std::string prop_name = name +
      //                  get_ion_name(ion); HDF5Tools::append_dataset< double
      //                  >(
      //                      group, prop_name, block_offset + offset,
      //                      scalar_props[scalar_index]);
      //                  ++scalar_index;
      //                }
      //              }
      //            } else if
      //            (DensityGridWriterFields::is_heating_property(property)) {
      //              for (int_fast32_t heating = 0; heating <
      //              NUMBER_OF_HEATINGTERMS;
      //                   ++heating) {
      //                if (fields.heatingterm_present(property, heating)) {
      //                  const std::string prop_name = name +
      //                  get_ion_name(heating); HDF5Tools::append_dataset<
      //                  double >(
      //                      group, prop_name, block_offset + offset,
      //                      scalar_props[scalar_index]);
      //                  ++scalar_index;
      //                }
      //              }
      //            } else {
      //              HDF5Tools::append_dataset< double >(group, name,
      //                                                  block_offset + offset,
      //                                                  scalar_props[scalar_index]);
      //              ++scalar_index;
      //            }
      //          }
      //        }
      //      }
      // WRITE ARRAYS TO FILE
    }
    block_offset += (*gridit).get_number_of_cells();
  }
  HDF5Tools::close_group(group);

  // close file
  HDF5Tools::close_file(file);
}
