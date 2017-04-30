//This file is part of Bertini 2.
//
//fixed_precision_tracker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fixed_precision_tracker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fixed_precision_tracker.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame


/**
\file fixed_precision_tracker.hpp

*/

#ifndef BERTINI_FIXED_PRECISION_TRACKER_HPP
#define BERTINI_FIXED_PRECISION_TRACKER_HPP

#include "bertini2/trackers/base_tracker.hpp"


namespace bertini{

	namespace tracking{

		using std::max;
		using std::min;
		using std::pow;

		using bertini::max;


	




		/** 
		\class FixedPrecisionTracker<prec>

		\brief Functor-like class for tracking paths on a system
		*/
		template<class DerivedT>
		class FixedPrecisionTracker : public Tracker<FixedPrecisionTracker<DerivedT>>
		{	
		public:

			using BaseComplexType = typename TrackerTraits<DerivedT>::BaseComplexType;
			using BaseRealType = typename TrackerTraits<DerivedT>::BaseRealType;

			using CT = BaseComplexType;
			using RT = BaseRealType;

			virtual ~FixedPrecisionTracker() = default;

			using EmitterType = FixedPrecisionTracker<DerivedT>;
			using Base = Tracker<FixedPrecisionTracker<DerivedT>>;

			using Config =  typename Base::Config;
			FORWARD_GET_CONFIGURED
			using Stepping = typename Base::Stepping;
			using Newton = typename Base::Newton;

			FixedPrecisionTracker(System const& sys);

			/**
			\brief An additional no-op call, provided for conformity of interface with AMP tracker in generic code.
			*/
			void PrecisionSetup(config::FixedPrecisionConfig<RT> const&)
			{ }


			Vec<CT> CurrentPoint() const override
			{
				return std::get<Vec<CT>>(this->current_space_);
			}

			virtual
			const RT LatestConditionNumber() const override
			{
				return std::get<RT>(this->condition_number_estimate_);
			}

			virtual
			const RT LatestErrorEstimate() const override
			{
				return std::get<RT>(this->error_estimate_);
			}

			virtual
			const RT LatestNormOfStep() const override
			{
				return std::get<RT>(this->norm_delta_z_);
			}


			void ResetCounters() const override
			{
				Base::ResetCountersBase();

				this->num_successful_steps_since_stepsize_increase_ = 0;
				// initialize to the frequency so guaranteed to compute it the first try 	
				this->num_steps_since_last_condition_number_computation_ = this->Get<Stepping>().frequency_of_CN_estimation;
			}


			/**
			\brief Ensure that number of steps, stepsize, and precision still ok.

			\return Success if ok to keep going, and a different code otherwise. 
			*/
			SuccessCode PreIterationCheck() const override
			{
				if (this->num_successful_steps_taken_ >= Get<Stepping>().max_num_steps)
					return SuccessCode::MaxNumStepsTaken;
				if (this->current_stepsize_ < Get<Stepping>().min_step_size)
					return SuccessCode::MinStepSizeReached;

				return SuccessCode::Success;
			}






			
			


			void PostTrackCleanup() const override
			{
				this->NotifyObservers(TrackingEnded<EmitterType>(*this));
			}

			/**
			\brief Copy from the internally stored current solution into a final solution.
			
			If preservation of precision is on, this function first returns to the initial precision.

			\param[out] solution_at_endtime The solution at the end time
			*/
			void CopyFinalSolution(Vec<CT> & solution_at_endtime) const override
			{

				// the current precision is the precision of the output solution point.

				unsigned num_vars = this->GetSystem().NumVariables();
				solution_at_endtime.resize(num_vars);
				for (unsigned ii=0; ii<num_vars; ii++)
				{
					solution_at_endtime(ii) = std::get<Vec<CT> >(this->current_space_)(ii);
				}

			}




			/**
			\brief Run an iteration of the tracker loop.

			Predict and correct, adjusting precision and stepsize as necessary.

			\return Success if the step was successful, and a non-success code if something went wrong, such as a linear algebra failure or AMP Criterion violation.
			*/
			SuccessCode TrackerIteration() const override;


			/**
			Check whether the path is going to infinity.
			*/
			SuccessCode CheckGoingToInfinity() const override
			{
				return Base::template CheckGoingToInfinity<CT>();
			}

			



			/**
			\brief Commit the next stepsize, and adjust internals.
			*/
			SuccessCode UpdateStepsize() const
			{
				this->SetStepSize(this->next_stepsize_);
				return SuccessCode::Success;
			}





			///////////
			//
			//  overrides for counter adjustment after a TrackerIteration()
			//
			////////////////

			/**
			\brief Increment and reset counters after a successful TrackerIteration()
			*/
			void OnStepSuccess() const override
			{
				Base::IncrementBaseCountersSuccess();
				this->NotifyObservers(SuccessfulStep<EmitterType >(*this));
			}

			/**
			\brief Increment and reset counters after a failed TrackerIteration()
			*/
			void OnStepFail() const override
			{
				Base::IncrementBaseCountersFail();
				this->num_successful_steps_since_stepsize_increase_ = 0;
				this->NotifyObservers(FailedStep<EmitterType >(*this));
			}



			void OnInfiniteTruncation() const override
			{
				this->NotifyObservers(InfinitePathTruncation<EmitterType>(*this));
			}

			//////////////
			//
			//


			/**
			\brief Wrapper function for calling the correct predictor.
			
			This function computes the next predicted space value, and sets some internals based on the prediction, such as the norm of the Jacobian.

			The real type and complex type must be commensurate.

			\param[out] predicted_space The result of the prediction
			\param current_space The current space point.
			\param current_time The current time value.
			\param delta_t The time differential for this step.  Allowed to be complex.
			*/
			SuccessCode Predict(Vec<CT> & predicted_space, 
								Vec<CT> const& current_space, 
								CT const& current_time, CT const& delta_t) const
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<CT>::Real, 
			              				typename Eigen::NumTraits<CT>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				RT& norm_J = std::get<RT>(this->norm_J_);
				RT& norm_J_inverse = std::get<RT>(this->norm_J_inverse_);
				RT& size_proportion = std::get<RT>(this->size_proportion_);
				RT& error_estimate = std::get<RT>(this->error_estimate_);
				RT& condition_number_estimate = std::get<RT>(this->condition_number_estimate_);

				return this->predictor_->Predict(
			                predicted_space,
							this->tracked_system_,
							current_space, current_time,
							delta_t,
							condition_number_estimate,
							this->num_steps_since_last_condition_number_computation_,
							Get<Stepping>().frequency_of_CN_estimation,
							RT(this->tracking_tolerance_));
			}



			/**
			\brief Run Newton's method.

			Wrapper function for calling Correct and getting the error estimates etc directly into the tracker object.

			\param corrected_space[out] The spatial result of the correction loop.
			\param current_space The start point in space for running the corrector loop.
			\param current_time The current time value.

			\return A SuccessCode indicating whether the loop was successful in converging in the max number of allowable newton steps, to the current path tolerance.
			*/
			SuccessCode Correct(Vec<CT> & corrected_space, 
								Vec<CT> const& current_space, 
								CT const& current_time) const
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RT>::Real, 
			              				typename Eigen::NumTraits<CT>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");


				RT& norm_J = std::get<RT>(this->norm_J_);
				RT& norm_J_inverse = std::get<RT>(this->norm_J_inverse_);
				RT& norm_delta_z = std::get<RT>(this->norm_delta_z_);
				RT& condition_number_estimate = std::get<RT>(this->condition_number_estimate_);


				return this->corrector_->Correct(corrected_space,
												this->tracked_system_,
												current_space,
												current_time, 
												RT(this->tracking_tolerance_),
												Get<Newton>().min_num_newton_iterations,
												Get<Newton>().max_num_newton_iterations);
			}



			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.  The tolerance is the tracking tolerance specified during Setup(...).


			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.

			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			SuccessCode RefineImpl(Vec<CT> & new_space,
								Vec<CT> const& start_point, CT const& current_time) const
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RT>::Real, 
			              				typename Eigen::NumTraits<CT>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				return this->corrector_->Correct(new_space,
							   this->tracked_system_,
							   start_point,
							   current_time, 
							   this->tracking_tolerance_,
							   Get<Newton>().min_num_newton_iterations,
							   Get<Newton>().max_num_newton_iterations);
			}






			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.


			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.
			\param tolerance The tolerance for convergence.  This is a tolerance on \f$\Delta x\f$, not on function residuals.
			\param max_iterations The maximum number of permitted Newton iterations.  

			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			SuccessCode RefineImpl(Vec<CT> & new_space,
								Vec<CT> const& start_point, CT const& current_time,
								RT const& tolerance, unsigned max_iterations) const
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RT>::Real, 
			              				typename Eigen::NumTraits<CT>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				return this->corrector_->Correct(new_space,
							   this->tracked_system_,
							   start_point,
							   current_time, 
							   tolerance,
							   1,
							   max_iterations);
			}


			/////////////////////////////////////////////
			//////////////////////////////////////
			/////////////////////////////
			////////////////////  data members stored in this class
			////////////
			//////
			//

			// no additional state variables needed for the FixedPrecision base tracker types

		};




		class DoublePrecisionTracker : public FixedPrecisionTracker<DoublePrecisionTracker>
		{
		public:
			using BaseComplexType = dbl;
			using BaseRealType = double;

			using EmitterType = typename TrackerTraits<DoublePrecisionTracker>::EventEmitterType;

			BERTINI_DEFAULT_VISITABLE()


			/**
			\brief Construct a tracker, associating to it a System.
			*/
			DoublePrecisionTracker(class System const& sys);


			DoublePrecisionTracker() = delete;

			virtual ~DoublePrecisionTracker() = default;


			unsigned CurrentPrecision() const override;


			/**
			\brief Set up the internals of the tracker for a fresh start.  

			Copies the start time, current stepsize, and start point.  Adjusts the current precision to match the precision of the start point.  Zeros counters.

			\param start_time The time at which to start tracking.
			\param end_time The time to which to track.
			\param start_point The space values from which to start tracking.
			*/
			SuccessCode TrackerLoopInitialization(BaseComplexType const& start_time,
			                               BaseComplexType const& end_time,
										   Vec<BaseComplexType> const& start_point) const override;


		private:

		}; // re: DoublePrecisionTracker


		class MultiplePrecisionTracker : public FixedPrecisionTracker<MultiplePrecisionTracker>
		{
		public:
			using BaseComplexType = mpfr;
			using BaseRealType = mpfr_float;

			using EmitterType = FixedPrecisionTracker<MultiplePrecisionTracker>;

			BERTINI_DEFAULT_VISITABLE()


			/**
			\brief Construct a tracker, associating to it a System.

			The precision of the tracker will be whatever the current default is.  The tracker cannot change its precision, and will require the default precision to be this precision whenever tracking is started.  That is, the precision is fixed.
			*/
			MultiplePrecisionTracker(class System const& sys);

			
			MultiplePrecisionTracker() = delete;

			virtual ~MultiplePrecisionTracker() = default;


			unsigned CurrentPrecision() const override;


			/**
			\brief Set up the internals of the tracker for a fresh start.  

			Copies the start time, current stepsize, and start point.  Adjusts the current precision to match the precision of the start point.  Zeros counters.

			\param start_time The time at which to start tracking.
			\param end_time The time to which to track.
			\param start_point The space values from which to start tracking.
			*/
			SuccessCode TrackerLoopInitialization(BaseComplexType const& start_time,
			                               BaseComplexType const& end_time,
										   Vec<BaseComplexType> const& start_point) const override;

			bool PrecisionSanityCheck() const;

		private:

			unsigned precision_;
		}; // re: MultiplePrecisionTracker


		extern template class FixedPrecisionTracker<DoublePrecisionTracker>;
		extern template class FixedPrecisionTracker<MultiplePrecisionTracker>;
	} // namespace tracking
} // namespace bertini


#endif



