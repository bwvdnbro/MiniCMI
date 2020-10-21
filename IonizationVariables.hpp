/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file IonizationVariables.hpp
 *
 * @brief Variables used in the ionization calculation.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef IONIZATIONVARIABLES_HPP
#define IONIZATIONVARIABLES_HPP

/**
 * @brief Variables used in the ionization calculation.
 */
class IonizationVariables {
private:
  /*! @brief Number density (in m^-3). */
  double _number_density;

  /*! @brief Temperature (in K). */
  double _temperature;

public:
  /**
   * @brief (Empty) constructor.
   */
  inline IonizationVariables()
      : _number_density(0.), _temperature(0.) {}

  /**
   * @brief Copy the contents of the given IonizationVariables instance into
   * this one.
   *
   * @param other Other IonizationVariables instance.
   */
  inline void copy_all(const IonizationVariables &other) {

    // single variables
    _number_density = other._number_density;
    _temperature = other._temperature;
  }

  /**
   * @brief Get the number density.
   *
   * @return Number density (in m^-3).
   */
  inline double get_number_density() const { return _number_density; }

  /**
   * @brief Set the number density.
   *
   * @param number_density New number density (in m^-3).
   */
  inline void set_number_density(const double number_density) {
    _number_density = number_density;
  }

  /**
   * @brief Get the temperature.
   *
   * @return Temperature (in K).
   */
  inline double get_temperature() const { return _temperature; }

  /**
   * @brief Set the temperature.
   *
   * @param temperature New temperature (in K).
   */
  inline void set_temperature(const double temperature) {
    _temperature = temperature;
  }
};

#endif // IONIZATIONVARIABLES_HPP
