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
 * @file HydroVariables.hpp
 *
 * @brief Hydro related variables.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef HYDROVARIABLES_HPP
#define HYDROVARIABLES_HPP

#include "CoordinateVector.hpp"

/**
 * @brief Hydro related variables.
 */
class HydroVariables {
private:
  /*! @brief Primitive variables. */
  double _primitives[5];

  /*! @brief Conserved variables. */
  double _conserved[5];

  /*! @brief Conserved variable changes. */
  double _delta_conserved[5];

  /*! @brief Gradients for the primitive variables. */
  CoordinateVector<> _primitive_gradients[5];

  /*! @brief Gravitational acceleration (in m s^-2). */
  CoordinateVector<> _gravitational_acceleration;

  /*! @brief External cooling and/or heating (in J s^-1). */
  double _energy_rate_term;

  /*! @brief External instanteneous energy terms (in J). */
  double _energy_term;

  /// TEST VARIABLES: these are used as an example

  /*! @brief Copy of the number density variable (in m^-3). */
  double _test_density;

  /*! @brief Sum of the _test_density values in all neighbouring cells
   *  (in m^-3). */
  double _test_neighbour_density_sum;

public:
  /**
   * @brief (Empty) constructor.
   */
  inline HydroVariables()
      : _primitives{0., 0., 0., 0., 0.}, _conserved{0., 0., 0., 0., 0.},
        _delta_conserved{0., 0., 0., 0., 0.}, _energy_rate_term(0.),
        _energy_term(0.),
        /* INITIALIZE TEST VARIABLES */
        _test_density(0.), _test_neighbour_density_sum(0.) {}

  /**
   * @brief Get read only access to the given component of the primitive
   * variables.
   *
   * @param index Index (0-4).
   * @return Read only access to the corresponding primitive variable component.
   */
  inline const double &primitives(uint_fast8_t index) const {
    return _primitives[index];
  }

  /**
   * @brief Get read/write access to the given component of the primitive
   * variables.
   *
   * @param index Index (0-4).
   * @return Read/write access to the corresponding primitive variable
   * component.
   */
  inline double &primitives(uint_fast8_t index) { return _primitives[index]; }

  /**
   * @brief Get the fluid density.
   *
   * @return Fluid density (in kg m^-3).
   */
  inline double get_primitives_density() const { return _primitives[0]; }

  /**
   * @brief Get the fluid velocity.
   *
   * @return Fluid velocity (in m s^-1).
   */
  inline CoordinateVector<> get_primitives_velocity() const {
    return CoordinateVector<>(_primitives[1], _primitives[2], _primitives[3]);
  }

  /**
   * @brief Get the fluid pressure.
   *
   * @return Fluid pressure (in kg m^-1 s^-2).
   */
  inline double get_primitives_pressure() const { return _primitives[4]; }

  /**
   * @brief Set the fluid density.
   *
   * @param density New fluid density (in kg m^-3).
   */
  inline void set_primitives_density(const double density) {
    _primitives[0] = density;
  }

  /**
   * @brief Set the fluid velocity.
   *
   * @param velocity New fluid velocity (in m s^-1).
   */
  inline void set_primitives_velocity(const CoordinateVector<> velocity) {
    _primitives[1] = velocity.x();
    _primitives[2] = velocity.y();
    _primitives[3] = velocity.z();
  }

  /**
   * @brief Set the fluid pressure.
   *
   * @param pressure New fluid pressure (in kg m^-1 s^-2).
   */
  inline void set_primitives_pressure(const double pressure) {
    _primitives[4] = pressure;
  }

  /**
   * @brief Get read only access to the given component of the conserved
   * variables.
   *
   * @param index Index (0-4).
   * @return Read only access to the corresponding conserved variable component.
   */
  inline const double &conserved(uint_fast8_t index) const {
    return _conserved[index];
  }

  /**
   * @brief Get read/write access to the given component of the conserved
   * variables.
   *
   * @param index Index (0-4).
   * @return Read/write access to the corresponding conserved variable
   * component.
   */
  inline double &conserved(uint_fast8_t index) { return _conserved[index]; }

  /**
   * @brief Get the fluid mass.
   *
   * @return Fluid mass (in kg).
   */
  inline double get_conserved_mass() const { return _conserved[0]; }

  /**
   * @brief Get the fluid momentum.
   *
   * @return Fluid momentum (in kg m s^-1).
   */
  inline CoordinateVector<> get_conserved_momentum() const {
    return CoordinateVector<>(_conserved[1], _conserved[2], _conserved[3]);
  }

  /**
   * @brief Get the fluid total energy.
   *
   * @return Fluid total energy (in kg m^2 s^-2).
   */
  inline double get_conserved_total_energy() const { return _conserved[4]; }

  /**
   * @brief Set the fluid mass.
   *
   * @param mass New fluid mass (in kg).
   */
  inline void set_conserved_mass(const double mass) { _conserved[0] = mass; }

  /**
   * @brief Set the fluid momentum.
   *
   * @param momentum New fluid momentum (in kg m s^-1).
   */
  inline void set_conserved_momentum(const CoordinateVector<> momentum) {
    _conserved[1] = momentum.x();
    _conserved[2] = momentum.y();
    _conserved[3] = momentum.z();
  }

  /**
   * @brief Set the fluid total energy.
   *
   * @param total_energy New fluid total energy (in kg m^2 s^-2).
   */
  inline void set_conserved_total_energy(const double total_energy) {
    _conserved[4] = total_energy;
  }

  /**
   * @brief Get read only access to the given component of the conserved
   * variable differences.
   *
   * @param index Index (0-4).
   * @return Read only access to the corresponding component of the conserved
   * variable differences.
   */
  inline const double &delta_conserved(uint_fast8_t index) const {
    return _delta_conserved[index];
  }

  /**
   * @brief Get read/write access to the given component of the conserved
   * variable differences.
   *
   * @param index Index (0-4).
   * @return Read/write access to the corresponding component of the conserved
   * variable differences.
   */
  inline double &delta_conserved(uint_fast8_t index) {
    return _delta_conserved[index];
  }

  /**
   * @brief Get read only access to the given component of the primitive
   * variable gradients.
   *
   * @param index Index (0-4).
   * @return Read only access to the corresponding component of the primitive
   * variable gradients.
   */
  inline const CoordinateVector<> &
  primitive_gradients(uint_fast8_t index) const {
    return _primitive_gradients[index];
  }

  /**
   * @brief Get read/write access to the given component of the primitive
   * variable gradients.
   *
   * @param index Index (0-4).
   * @return Read/write access to the corresponding component of the primitive
   * variable gradients.
   */
  inline CoordinateVector<> &primitive_gradients(uint_fast8_t index) {
    return _primitive_gradients[index];
  }

  /**
   * @brief Get the gravitational acceleration.
   *
   * @return Gravitational acceleration (in m s^-2).
   */
  inline const CoordinateVector<> get_gravitational_acceleration() const {
    return _gravitational_acceleration;
  }

  /**
   * @brief Set the gravitational acceleration.
   *
   * @param gravitational_acceleration Gravitational acceleration (in m s^-2).
   */
  inline void set_gravitational_acceleration(
      const CoordinateVector<> gravitational_acceleration) {
    _gravitational_acceleration = gravitational_acceleration;
  }

  /**
   * @brief Get the energy rate term.
   *
   * @return Energy rate term (in J s^-1).
   */
  inline double get_energy_rate_term() const { return _energy_rate_term; }

  /**
   * @brief Set the energy rate term.
   *
   * @param energy_rate_term Energy rate term (in J s^-1).
   */
  inline void set_energy_rate_term(const double energy_rate_term) {
    _energy_rate_term = energy_rate_term;
  }

  /**
   * @brief Get the energy term.
   *
   * @return Energy term (in J).
   */
  inline double get_energy_term() const { return _energy_term; }

  /**
   * @brief Set the energy term.
   *
   * @param energy_term Energy term (in J).
   */
  inline void set_energy_term(const double energy_term) {
    _energy_term = energy_term;
  }

  /**
   * @brief Copy the contents of the given HydroVariables instance into
   * this one.
   *
   * @param other Other HydroVariables instance.
   */
  inline void copy_all(const HydroVariables &other) {

    for (uint_fast8_t i = 0; i < 5; ++i) {
      _primitives[i] = other._primitives[i];
      _conserved[i] = other._conserved[i];
      _delta_conserved[i] = other._delta_conserved[i];
      _primitive_gradients[i] = other._primitive_gradients[i];
      _gravitational_acceleration[i] = other._gravitational_acceleration[i];
    }

    _energy_rate_term = other._energy_rate_term;
    _energy_term = other._energy_term;
  }

  /// TEST VARIABLE GETTERS/SETTERS

  inline double get_test_density() const { return _test_density; }

  inline void set_test_density(const double test_density) {
    _test_density = test_density;
  }

  inline double get_test_neighbour_density_sum() const {
    return _test_neighbour_density_sum;
  }

  inline void
  set_test_neighbour_density_sum(const double test_neighbour_density_sum) {
    _test_neighbour_density_sum = test_neighbour_density_sum;
  }
};

#endif // HYDROVARIABLES_HPP
