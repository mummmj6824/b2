//This file is part of Bertini 2.
//
//amp_tracker.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_tracker.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_tracker.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame


/**
\file amp_tracker.hpp

\brief 

\brief Contains the AMPTracker type, and necessary functions.
*/

#ifndef BERTINI_AMP_TRACKER_HPP
#define BERTINI_AMP_TRACKER_HPP

#include "bertini2/trackers/base_tracker.hpp"


namespace bertini{

	namespace tracking{

		
		using std::max;
		using std::min;
		using std::pow;

		using bertini::max;


		inline
		mpfr_float MinTimeForCurrentPrecision(unsigned precision)
		{
			if (precision==DoublePrecision())
				return max( mpfr_float(pow( mpfr_float(10), -long(precision)+3)), mpfr_float("1e-150"));
			else
				return pow( mpfr_float(10), -long(precision)+3);
		}

		inline
		mpfr_float MinStepSizeForPrecision(unsigned precision)
		{
			return MinTimeForCurrentPrecision(precision);
		}


		/**
		 \brief Compute the cost function for arithmetic versus precision.

		 From \cite AMP2, \f$C(P)\f$.  As currently implemented, this is 
		 \f$ 10.35 + 0.13 P \f$, where P is the precision.
		
		 This function tells you the relative cost of arithmetic at a given precision. 

		 \todo Recompute this cost function for boost::multiprecision::mpfr_float
		*/
		inline
		mpfr_float ArithmeticCost(unsigned precision)
		{
			if (precision==DoublePrecision())
				return 1;
			else
				return mpfr_float("10.35") + mpfr_float("0.13") * precision;
		}




		/**
		 \Compute a stepsize satisfying AMP Criterion B with a given precision
		*/
		inline
		mpfr_float StepsizeSatisfyingCriterionB(unsigned precision,
										mpfr_float const& criterion_B_rhs,
										unsigned num_newton_iterations,
										unsigned predictor_order = 0)
		{
			return pow(mpfr_float(10), -( criterion_B_rhs - precision )*num_newton_iterations/(predictor_order+1));
		}
		

		/**
		\brief The minimum number of digits based on a log10 of a stepsize.

		\param log_of_stepsize log10 of a stepsize
		\param current_time The current time
		*/
		inline
		unsigned MinDigitsForLogOfStepsize(mpfr_float const& log_of_stepsize, mpfr const& current_time)
		{
			return mpfr_float(ceil(log_of_stepsize) + ceil(log10(abs(current_time))) + 2).convert_to<int>();
		}

		inline
		unsigned MinDigitsForStepsizeInterval(mpfr_float const& min_stepsize, mpfr_float const& max_stepsize, mpfr const& current_time)
		{
			return max(MinDigitsForLogOfStepsize(-log10(min_stepsize),current_time), MinDigitsForLogOfStepsize(-log10(max_stepsize),current_time));
		}
		

		/**
		 \brief Compute precision and stepsize minimizing the ArithmeticCost() of tracking.
		 
		 For a given range of precisions, an old stepsize, and a maximum stepsize, the Cost of tracking is computed, and a minimizer found.  
		
		 This function is used in the AMPTracker tracking loop, both in case of successful steps and in Criterion errors.


		 \param[out] new_precision The minimizing precision.
		 \param[out] new_stepsize The minimizing stepsize.
		 \param[in] min_precision The minimum considered precision.
		 \param[in] old_stepsize The previously used stepsize.
		 \param[in] max_precision The maximum considered precision.
		 \param[in] max_stepsize The maximum permitted stepsize.  
		 \param[in] criterion_B_rhs The right hand side of CriterionB from \cite AMP1, \cite AMP2
		 \param num_newton_iterations The number of allowed Newton corrector iterations.
		 \param AMP_config The configuration of AMP settings for tracking.
		 \param predictor_order The order of the predictor being used.  This is the order itself, not the order of the error estimate.

		 \see ArithmeticCost
		*/
		void MinimizeTrackingCost(unsigned & new_precision, mpfr_float & new_stepsize, 
						  unsigned min_precision, mpfr_float const& old_stepsize,
						  unsigned max_precision, mpfr_float const& max_stepsize,
						  mpfr_float const& criterion_B_rhs,
						  unsigned num_newton_iterations,
						  unsigned predictor_order = 0);








		/** 
		\class AMPTracker

		\brief Functor-like class for tracking paths on a system
		
		
		## Explanation
		
		The bertini::AMPTracker class enables tracking using Adaptive Multiple Precision on an arbitrary square homotopy.  

		The intended usage is to:

		1. Create a system, and instantiate some settings.
		2. Create an AMPTracker, associating it to the system you are going to solve or track on.
		3. Run AMPTracker::Setup and AMPTracker::PrecisionSetup, getting the settings in line for tracking.
		4. Repeatedly, or as needed, call the AMPTracker::TrackPath function, feeding it a start point, and start and end times.  The initial precision is that of the start point.  

		Working precision and stepsize are adjusted automatically to get around nearby singularities to the path.  If the endpoint is singular, this may very well fail, as prediction and correction get more and more difficult with proximity to singularities.  

		The TrackPath method is intended to allow the user to track to nonsingular endpoints, or to an endgame boundary, from which an appropriate endgame will be called.
		
		## Some notes

		The AMPTracker has internal state.  That is, it stores the current state of the path being tracked as data members.  After tracking has concluded, these statistics may be extracted.  If you need additional accessors for these data, contact the software authors.

		The class additionally uses the Boost.Log library for logging.  At time of this writing, the trivial logger is being used.  As development continues we will move toward using a more sophisticated logging model.  Suggestions are welcome.
		
		This class, like the other Tracker classes, uses mutable members to store the current state of the tracker.  The tracking tolerances, settings, etc, remain constant throughout a track, but the internal state such as the current time or space values, will change.  You can also expect the precision of the tracker to differ after a certain calls, too.
		
		## Example Usage
		
		Below we demonstrate a basic usage of the AMPTracker class to track a single path.  

		The pattern is as described above: create an instance of the class, feeding it the system to be tracked, and some configuration.  Then, use the tracker to track paths of the system.

		\code{.cpp}
		DefaultPrecision(30); // set initial precision.  This is not strictly necessary.

		using namespace bertini::tracking;

		// 1. Create the system
		Var x = MakeVariable("x");
		Var y = MakeVariable("y");
		Var t = MakeVariable("t");

		System sys;

		VariableGroup v{x,y};

		sys.AddFunction(pow(x,2) + (1-t)*x - 1);
		sys.AddFunction(pow(y,2) + (1-t)*x*y - 2);
		sys.AddPathVariable(t);
		sys.AddVariableGroup(v);

		auto AMP = bertini::tracking::config::AMPConfigFrom(sys);
		
		//  2. Create the Tracker object, associating the system to it.
		bertini::tracking::AMPTracker tracker(sys);

		config::Stepping stepping_preferences;
		config::Newton newton_preferences;

		// 3. Get the settings into the tracker 
		tracker.Setup(config::Predictor::Euler,
		              	mpfr_float("1e-5"),
						mpfr_float("1e5"),
						stepping_preferences,
						newton_preferences);

		tracker.PrecisionSetup(AMP);
	
		//  4. Create a start and end time.  These are complex numbers.
		mpfr t_start("1.0");
		mpfr t_end("0");
		
		//  5. Create a start point, and container for the end point.
		Vec<mpfr> start_point(2);
		start_point << mpfr("1"), mpfr("1.414");  // set the value of the start point.  This is Eigen syntax.

		Vec<mpfr> end_point;

		// 6. actually do the tracking
		SuccessCode tracking_success = tracker.TrackPath(end_point,
		                  t_start, t_end, start_point);
		
		// 7. and then onto whatever processing you are doing to the computed point.
		\endcode
		
		If this documentation is insufficient, please contact the authors with suggestions, or get involved!  Pull requests welcomed.
		
		## Testing

		* Test suite driving this class: AMP_tracker_basics.
		* File: test/tracking_basics/tracker_test.cpp
		* Functionality tested: Can use an AMPTracker to track in various situations, including tracking on nonsingular paths.  Also track to a singularity on the square root function.  Furthermore test that tracking fails to start from a singular start point, and tracking fails if a singularity is directly on the path being tracked.

		*/
		class AMPTracker : public Tracker<AMPTracker>
		{
			friend class Tracker<AMPTracker>;
		public:
			BERTINI_DEFAULT_VISITABLE()
			
			typedef Tracker<AMPTracker> Base;
			typedef typename TrackerTraits<AMPTracker>::EventEmitterType EmitterType;

			enum UpsampleRefinementOption
			{
			   upsample_refine_off  = 0,
			   upsample_refine_on   = 1
			};
		

			/**
			\brief Construct an Adaptive Precision tracker, associating to it a System.
			*/
			AMPTracker(class System const& sys) : Tracker(sys), current_precision_(DefaultPrecision())
			{	
				Set<PrecConf>(config::AMPConfigFrom(sys));
			}
			

			/**
			\brief Special additional setup call for the AMPTracker, selecting the config for adaptive precision.
			*/
			void PrecisionSetup(config::AdaptiveMultiplePrecisionConfig const& AMP_config)
			{
				Set<PrecConf>(AMP_config);
			}


			const unsigned GetCurrentPrecision() const
			{
				return current_precision_;
			}

			
			/**
			\brief Switch preservation of precision after tracking on / off

			By default, precision is preserved after tracking, so the precision of the ambient workspace is returned to its previous state once Adaptive Precision tracking is done.
			*/
			void PrecisionPreservation(bool should_preseve_precision)
			{
				preserve_precision_ = should_preseve_precision;
			}

			
			virtual ~AMPTracker() = default;


			Vec<mpfr> CurrentPoint() const override
			{
				if (this->CurrentPrecision()==DoublePrecision())
				{
					const auto& curr_vector = std::get<Vec<dbl>>(this->current_space_);
					Vec<mpfr> returnme(NumVariables());
					for (unsigned ii = 0; ii < NumVariables(); ++ii)
					{
						returnme(ii) = mpfr(curr_vector(ii));
					}
					return returnme;
				}
				else
					return std::get<Vec<mpfr>>(this->current_space_);
			}

			virtual
			const mpfr_float LatestConditionNumber() const override
			{
				if (this->CurrentPrecision()!=DoublePrecision())
					return std::get<mpfr_float>(this->condition_number_estimate_);
				else
					return mpfr_float(std::get<double>(this->condition_number_estimate_));
			}

			virtual
			const mpfr_float LatestErrorEstimate() const override
			{
				if (this->CurrentPrecision()!=DoublePrecision())
					return std::get<mpfr_float>(this->error_estimate_);
				else
					return mpfr_float(std::get<double>(this->error_estimate_));
			}

			virtual
			const mpfr_float LatestNormOfStep() const override
			{
				if (this->CurrentPrecision()!=DoublePrecision())
					return std::get<mpfr_float>(this->norm_delta_z_);
				else
					return mpfr_float(std::get<double>(this->norm_delta_z_));
			}



		private:

			/**
			\brief Set up the internals of the tracker for a fresh start.  

			Copies the start time, current stepsize, and start point.  Adjusts the current precision to match the precision of the start point.  Zeros counters.

			\param start_time The time at which to start tracking.
			\param end_time The time to which to track.
			\param start_point The space values from which to start tracking.
			*/
			SuccessCode TrackerLoopInitialization(mpfr const& start_time,
			                               mpfr const& end_time,
										   Vec<mpfr> const& start_point) const override;


			// SuccessCode TrackerLoopInitialization(dbl const& start_time,
			//                                dbl const& end_time,
			// 							   Vec<dbl> const& start_point) const override
			// {
			// 	NotifyObservers(Initializing<AMPTracker,dbl>(*this,start_time, end_time, start_point));

			// 	// set up the master current time and the current step size
			// 	initial_precision_ = Precision(DoublePrecision());
			// 	DefaultPrecision(initial_precision_);
			// 	current_time_.precision(initial_precision_);
			// 	current_time_ = mpfr(start_time);

			// 	current_stepsize_.precision(initial_precision_);
			// 	if (reinitialize_stepsize_)
			// 		SetStepSize(min(double(Get<Stepping>().initial_step_size),abs(start_time-end_time)/Get<Stepping>().min_num_steps));

			// 	// populate the current space value with the start point, in appropriate precision
				
			// 	DoubleToDouble( start_point);
			// 	ChangePrecision<upsample_refine_off>(DoublePrecision());
			// 	ResetCounters();
			// 	return InitialRefinement();
			// }


			void ResetCounters() const override;

			/** 
			\brief Run an initial refinement of the start point, to ensure in high enough precision to start.

			\return Whether initial refinement was successful.  
			*/
			SuccessCode InitialRefinement() const;




			/**
			\brief Ensure that number of steps, stepsize, and precision still ok.

			\return Success if ok to keep going, and a different code otherwise. 
			*/
			SuccessCode PreIterationCheck() const override;



			
			


			void PostTrackCleanup() const override;

			/**
			\brief Copy from the internally stored current solution into a final solution.
			
			If preservation of precision is on, this function first returns to the initial precision.

			\param[out] solution_at_endtime The solution at the end time
			*/
			void CopyFinalSolution(Vec<mpfr> & solution_at_endtime) const override;


			// void CopyFinalSolution(Vec<dbl> & solution_at_endtime) const override
			// {
			// 	if (preserve_precision_)
			// 		ChangePrecision(initial_precision_);

			// 	// the current precision is the precision of the output solution point.
			// 	if (current_precision_==DoublePrecision())
			// 	{
			// 		solution_at_endtime = std::get<Vec<dbl> >(current_space_);
			// 	}
			// 	else
			// 	{
			// 		unsigned num_vars = GetSystem().NumVariables();
			// 		solution_at_endtime.resize(num_vars);
			// 		for (unsigned ii=0; ii<num_vars; ii++)
			// 			solution_at_endtime(ii) = dbl(std::get<Vec<mpfr> >(current_space_)(ii));

			// 	}
			// }

			/**
			\brief Run an iteration of the tracker loop.

			Predict and correct, adjusting precision and stepsize as necessary.

			\return Success if the step was successful, and a non-success code if something went wrong, such as a linear algebra failure or AMP Criterion violation.
			*/
			SuccessCode TrackerIteration() const override;


			/**
			\brief Run an iteration of AMP tracking.

			\return SuccessCode indicating whether the iteration was successful.
			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode TrackerIteration() const;


			/**
			Check whether the path is going to infinity.
			*/
			SuccessCode CheckGoingToInfinity() const override
			{
				if (current_precision_ == DoublePrecision())
					return Base::CheckGoingToInfinity<dbl>();
				else
					return Base::CheckGoingToInfinity<mpfr>();
			}

			/**
			\brief Commit the next precision and stepsize, and adjust internals.

			\tparam refine_if_necessary Flag indicating whether to refine if precision if increasing.

			\see UpsampleRefinementOption
			*/
			template <UpsampleRefinementOption refine_if_necessary = upsample_refine_off>
			SuccessCode UpdatePrecisionAndStepsize() const
			{
				SetStepSize(next_stepsize_);
				return ChangePrecision<refine_if_necessary>(next_precision_);
			}







			/**
			\brief Increase stepsize or decrease precision, because of consecutive successful steps.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			If the most recent step was successful, maybe adjust down precision and up stepsize.  

			The number of consecutive successful steps is recorded as state in this class, and if this number exceeds a user-determine threshold, the precision or stepsize are allowed to favorably change.  If not, then precision can only go up or remain the same.  Stepsize can only decrease.  These changes depend on the AMP criteria and current tracking tolerance.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode AdjustAMPStepSuccess() const;







			/**
			\brief The convergence_error function from \cite AMP2.  
	
			Adjust precision and stepsize due to Newton's method failure to converge in the maximum number of steps.

			Decrease stepsize, then increase precision until the new stepsize is greater than a precision-dependent minimum stepsize.  
			
			This function is used in the AMPTracker loop, when Newton's method fails to converge.

			\see MinStepSizeForPrecision
			*/
			void NewtonConvergenceError() const;







			/**
			\brief Adjust step size and precision due to AMP Criterion violation
			
			This function adjusts internals to the tracker object.
			
			Precision is REQUIRED to increase at least one increment.  Stepsize is REQUIRED to decrease.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template<typename ComplexType, typename RealType>
			void AMPCriterionError() const;





			/**
			\brief Get the raw right-hand side of Criterion B based on current state.
			*/
			template<typename ComplexType, typename RealType>
			RealType B_RHS() const
			{	
				return max(amp::CriterionBRHS(std::get<RealType>(norm_J_), 
				           					  std::get<RealType>(norm_J_inverse_), 
				           					  Get<Newton>().max_num_newton_iterations, 
				           					  RealType(tracking_tolerance_), 
				           					  std::get<RealType>(size_proportion_), 
				           					  Get<PrecConf>()),
				            RealType(0));
			}




			/**
			\brief Get the right hand side of Criterion B based on current state.
			
			Returns the larger of the right hand side of B, or 1.

			Uses:

			* the most recent estimate on the norm of \f$J\f$,
			* the most recent estimate of the norm of \f$J^{-1}\f$,
			* the number of allowed newton iterations,
			* the tracking tolerance,
			* the size_proportion of the latest prediction,
			* the AMP configuration.

			*/
			template<typename ComplexType, typename RealType>
			unsigned DigitsB() const
			{	
				return unsigned(B_RHS<ComplexType, RealType>());
			}



			

			/**
			\brief Get the raw right-hand side of Criterion B based on current state.
			*/
			template<typename ComplexType, typename RealType>
			RealType C_RHS() const
			{	
				return max(amp::CriterionCRHS(std::get<RealType>(norm_J_inverse_), 
				                              std::get<Vec<ComplexType> > (current_space_).norm(), 
				                              RealType(tracking_tolerance_), 
				                              Get<PrecConf>()),
				           RealType(0));
			}



			/**
			\brief Get the right hand side of Criterion C based on current state.
			
			Returns the larger of the right hand side of C, or 1.

			Uses:

			* the most recent estimate of the norm of \f$J^{-1}\f$,
			* the most recent current space point,
			* the current tracking tolerance,
			* the AMP configuration.

			*/
			template<typename ComplexType, typename RealType>
			unsigned DigitsC() const
			{	
				return unsigned(C_RHS<ComplexType, RealType>());
			}





			/**
			\brief Get the minimum required precision based on current state.

			The current state determining the minimum required precision is:

			* the norm of the Jacobian,
			* the norm of the inverse of the Jacobian,
			* the size_proportion, related to stepsize and predictor order,
			* the norm of the current solution. 
			* the tracking tolerance.

			The min digits needed is the maximum of 

			* Criterion B, 
			* Criterion C, and 
			* the digits required by the tracking tolerance.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template <typename ComplexType, typename RealType>
			unsigned MinRequiredPrecision_BCTol() const
			{
				return max(DigitsB<ComplexType,RealType>(), 
				           DigitsC<ComplexType,RealType>(), 
				           digits_tracking_tolerance_, 
				           DoublePrecision()); 
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
				Tracker::IncrementBaseCountersSuccess();
				NotifyObservers(SuccessfulStep<EmitterType>(*this));
			}

			/**
			\brief Increment and reset counters after a failed TrackerIteration()
			*/
			void OnStepFail() const override
			{
				Tracker::IncrementBaseCountersFail();
				num_successful_steps_since_precision_decrease_ = 0;
				num_successful_steps_since_stepsize_increase_ = 0;
				NotifyObservers(FailedStep<EmitterType>(*this));
			}



			void OnInfiniteTruncation() const override
			{
				NotifyObservers(InfinitePathTruncation<EmitterType>(*this));
			}


			////////////////
			//
			//       Predict and Correct functions
			//
			/////////////////


			


/**
			\brief Wrapper function for calling the correct predictor.
			
			This function computes the next predicted space value, and sets some internals based on the prediction, such as the norm of the Jacobian.

			The real type and complex type must be commensurate.

			\param[out] predicted_space The result of the prediction
			\param current_space The current space point.
			\param current_time The current time value.
			\param delta_t The time differential for this step.  Allowed to be complex.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.
			*/
			template<typename ComplexType, typename RealType, typename Derived>
			SuccessCode Predict(Vec<ComplexType> & predicted_space, 
								const Eigen::MatrixBase<Derived>& current_space,
								ComplexType const& current_time, ComplexType const& delta_t) const
			{
				
								
				
				
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");
				static_assert(std::is_same<typename Derived::Scalar, ComplexType>::value, "scalar types must match");

				RealType& norm_J = std::get<RealType>(norm_J_);
				RealType& norm_J_inverse = std::get<RealType>(norm_J_inverse_);
				RealType& size_proportion = std::get<RealType>(size_proportion_);
				RealType& error_estimate = std::get<RealType>(error_estimate_);
				RealType& condition_number_estimate = std::get<RealType>(condition_number_estimate_);

				if (predictor_->HasErrorEstimate())
					return predictor_->Predict<ComplexType,RealType>(predicted_space,
									error_estimate,
									size_proportion,
									norm_J,
									norm_J_inverse,
									tracked_system_,
									current_space, current_time, 
									delta_t,
									condition_number_estimate,
									num_steps_since_last_condition_number_computation_, 
									Get<Stepping>().frequency_of_CN_estimation, 
									RealType(tracking_tolerance_),
									Get<PrecConf>());
				else
					return predictor_->Predict<ComplexType,RealType>(predicted_space,
									size_proportion,
									norm_J,
									norm_J_inverse,
									tracked_system_,
									current_space, current_time, 
									delta_t,
									condition_number_estimate,
									num_steps_since_last_condition_number_computation_, 
									Get<Stepping>().frequency_of_CN_estimation, 
									RealType(tracking_tolerance_),
									Get<PrecConf>());
			}



			/**
			\brief Run Newton's method.

			Wrapper function for calling Correct and getting the error estimates etc directly into the tracker object.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			\param corrected_space[out] The spatial result of the correction loop.
			\param current_space The start point in space for running the corrector loop.
			\param current_time The current time value.

			\return A SuccessCode indicating whether the loop was successful in converging in the max number of allowable newton steps, to the current path tolerance.
			*/
			template<typename ComplexType, typename RealType>
			SuccessCode Correct(Vec<ComplexType> & corrected_space, 
								Vec<ComplexType> const& current_space, 
								ComplexType const& current_time) const
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");


				RealType& norm_J = std::get<RealType>(norm_J_);
				RealType& norm_J_inverse = std::get<RealType>(norm_J_inverse_);
				RealType& norm_delta_z = std::get<RealType>(norm_delta_z_);
				RealType& condition_number_estimate = std::get<RealType>(condition_number_estimate_);


				return corrector_->Correct(corrected_space,
									norm_delta_z,
									norm_J,
									norm_J_inverse,
									condition_number_estimate,
									tracked_system_,
									current_space,
									current_time,
									RealType(tracking_tolerance_),
									Get<Newton>().min_num_newton_iterations,
									Get<Newton>().max_num_newton_iterations,
									Get<PrecConf>());
			}



			/**
			\brief Run Newton's method from the currently stored time and space value, in current precision.

			This overwrites the value of the current space, if successful, with the refined value.

			\return Whether the refinement was successful.
			*/
			SuccessCode RefineStoredPoint() const
			{
				SuccessCode code;
				if (current_precision_==DoublePrecision())
				{
					code = RefineImpl<dbl>(std::get<Vec<dbl> >(temporary_space_),std::get<Vec<dbl> >(current_space_), dbl(current_time_));
					if (code == SuccessCode::Success)
						std::get<Vec<dbl> >(current_space_) = std::get<Vec<dbl> >(temporary_space_);
				}
				else
				{
					code = RefineImpl<mpfr>(std::get<Vec<mpfr> >(temporary_space_),std::get<Vec<mpfr> >(current_space_), current_time_);
					if (code == SuccessCode::Success)
						std::get<Vec<mpfr> >(current_space_) = std::get<Vec<mpfr> >(temporary_space_);
				}
				return code;
			}



			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.  The tolerance is the tracking tolerance specified during Setup(...).

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.

			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			template <typename ComplexType>
			SuccessCode RefineImpl(Vec<ComplexType> & new_space,
								Vec<ComplexType> const& start_point, ComplexType const& current_time) const
			{
				using RealType = typename Eigen::NumTraits<ComplexType>::Real;
				
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");

				auto target_precision = Precision(current_time);
				assert(Precision(start_point)==target_precision);
				ChangePrecision(target_precision);
				Precision(new_space,target_precision);
				
				RealType& norm_J = std::get<RealType>(norm_J_);
				RealType& norm_J_inverse = std::get<RealType>(norm_J_inverse_);
				RealType& norm_delta_z = std::get<RealType>(norm_delta_z_);
				RealType& condition_number_estimate = std::get<RealType>(condition_number_estimate_);


				return corrector_->Correct(new_space,
										   norm_delta_z,
										   norm_J,
										   norm_J_inverse,
										   condition_number_estimate,
										   tracked_system_,
										   start_point,
										   current_time,
										   RealType(tracking_tolerance_),
										   Get<Newton>().min_num_newton_iterations,
										   Get<Newton>().max_num_newton_iterations,
										   Get<PrecConf>());
			}






			/**
			\brief Run Newton's method from a start point with a current time.  

			Returns new space point by reference, as new_space.  Operates at current precision.

			\tparam ComplexType The complex number type.
			\tparam RealType The real number type.

			\param[out] new_space The result of running the refinement.
			\param start_point The base point for running Newton's method.
			\param current_time The current time value.
			\param tolerance The tolerance for convergence.  This is a tolerance on \f$\Delta x\f$, not on function residuals.
			\param max_iterations The maximum allowable number of iterations to perform.
			\return Code indicating whether was successful or not.  Regardless, the value of new_space is overwritten with the correction result.
			*/
			template <typename ComplexType, typename RealType>
			SuccessCode RefineImpl(Vec<ComplexType> & new_space,
								Vec<ComplexType> const& start_point, ComplexType const& current_time,
								RealType const& tolerance, unsigned max_iterations) const
			{
				static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
			              				typename Eigen::NumTraits<ComplexType>::Real>::value,
			              				"underlying complex type and the type for comparisons must match");
				auto target_precision = Precision(current_time);
				assert(Precision(start_point)==target_precision);
				ChangePrecision(target_precision);
				Precision(new_space,target_precision);

				using R = typename Eigen::NumTraits<ComplexType>::Real;

				R& norm_J = std::get<R>(norm_J_);
				R& norm_J_inverse = std::get<R>(norm_J_inverse_);
				R& norm_delta_z = std::get<R>(norm_delta_z_);
				R& condition_number_estimate = std::get<R>(condition_number_estimate_);

				return corrector_->Correct(new_space,
							   norm_delta_z,
								norm_J,
								norm_J_inverse,
								condition_number_estimate,
							   tracked_system_,
							   start_point,
							   current_time, 
							   tolerance,
							   1,
							   max_iterations,
							   Get<PrecConf>());
			}














			/////////////////
			//
			//  Functions for converting between precision types
			//
			///////////////////////

		public:
			/**
			Change precision of tracker to next_precision.  Converts the internal temporaries, and adjusts precision of system. Then refines if necessary.

			If the new precision is higher than current precision, a refine step will be called, which runs Newton's method.  This may fail, leaving the tracker in a state with higher precision internals, but garbage digits after the previously known digits.

			\param new_precision The precision to change to.
			\return SuccessCode indicating whether the change was successful.  If the precision increases, and the refinement loop fails, this could be not Success.  Changing down is guaranteed to succeed.
			*/
			template <UpsampleRefinementOption refine_if_necessary = upsample_refine_off>
			SuccessCode ChangePrecision(unsigned new_precision) const
			{
				if (new_precision==current_precision_) // no op
					return SuccessCode::Success;

				NotifyObservers(PrecisionChanged<EmitterType>(*this,current_precision_,new_precision));
				

				bool upsampling_needed = new_precision > current_precision_;
				// reset the counter for estimating the condition number.  
				num_steps_since_last_condition_number_computation_ = this->Get<Stepping>().frequency_of_CN_estimation;

				if (new_precision==DoublePrecision() && current_precision_>DoublePrecision())
				{
					// convert from multiple precision to double precision
					MultipleToDouble();
				}
				else if(new_precision > DoublePrecision() && current_precision_ == DoublePrecision())
				{
					// convert from double to multiple precision
					DoubleToMultiple(new_precision);
				}
				else
				{
					MultipleToMultiple(new_precision);
				}



				#ifndef BERTINI_DISABLE_ASSERTS
				assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
				#endif


				if (refine_if_necessary && upsampling_needed)
					return RefineStoredPoint();
				else
					return SuccessCode::Success;
			}





			private:
			
			/**
			\brief Converts from double to double

			Copies a multiple-precision into the double storage vector, and changes precision of the time and delta_t.

			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void DoubleToDouble(Vec<dbl> const& source_point) const;

			/**
			\brief Converts from multiple to double

			Changes the precision of the internal temporaries to double precision
			*/
			void DoubleToDouble() const
			{
				DoubleToDouble(std::get<Vec<dbl> >(current_space_));
			}
			


			/**
			\brief Converts from multiple to double

			Copies a multiple-precision into the double storage vector, and changes precision of the time and delta_t.

			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void MultipleToDouble(Vec<mpfr> const& source_point) const;
			/**
			\brief Converts from multiple to double

			Changes the precision of the internal temporaries to double precision
			*/
			void MultipleToDouble() const
			{
				MultipleToDouble(std::get<Vec<mpfr> >(current_space_));
			}




			/**
			\brief Converts from double to multiple

			Copies a double-precision into the multiple-precision storage vector.  

			You should call Refine after this to populate the new digits with non-garbage data.
	
			\param new_precision The new precision.
			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void DoubleToMultiple(unsigned new_precision, Vec<dbl> const& source_point) const;



			/**
			\brief Converts from double to multiple

			Copies the double-precision temporaries into the multiple-precision temporaries.  You should call Newton after this to populate the new digits with non-garbage data.

			\param new_precision The new precision.
			*/
			void DoubleToMultiple(unsigned new_precision) const
			{
				DoubleToMultiple( new_precision, std::get<Vec<dbl> >(current_space_));
			}





			/**
			\brief Converts from multiple to different precision multiple precision

			Copies a multiple-precision into the multiple-precision storage vector, and changes precision of the time and delta_t.
			Also resets counter so have to re-compute the condition number on next step attempt.

			\param new_precision The new precision.
			\param source_point The point into which to copy to the internally stored current space point.
			*/
			void MultipleToMultiple(unsigned new_precision, Vec<mpfr> const& source_point) const;


			/**
			\brief Converts from multiple to different precision multiple precision

			Changes the precision of the internal temporaries to desired precision

			\param new_precision The new precision.
			*/
			void MultipleToMultiple(unsigned new_precision) const
			{
				MultipleToMultiple( new_precision, std::get<Vec<mpfr> >(current_space_));
			}


			/**
			\brief Change precision of all temporary internal state variables.

			This excludes those which cannot be re-written without copying -- the current space point most notably.

			\brief new_precision The new precision to adjust to.
			*/
			void AdjustTemporariesPrecision(unsigned new_precision) const;



			
			/**
			\brief Ensure that all internal state is in uniform precision.

			\return True if all internal state variables are in the same precision as current_precision_, false otherwise.
			*/
			bool PrecisionSanityCheck() const;




			/////////////////////////////////////////////
			//////////////////////////////////////
			/////////////////////////////
			////////////////////  data members stored in this class
			////////////
			//////
			//





			////////////
			// state variables
			/////////////
			bool preserve_precision_ = false; ///< Whether the tracker should change back to the initial precision after tracking paths.

			mutable unsigned previous_precision_; ///< The previous precision of the tracker.
			mutable unsigned current_precision_; ///< The current precision of the tracker, the system, and all temporaries.
			mutable unsigned next_precision_; ///< The next precision
			mutable unsigned num_precision_decreases_; ///< The number of times precision has decreased this track.
			mutable unsigned initial_precision_; ///< The precision at the start of tracking.
			mutable unsigned num_successful_steps_since_precision_decrease_; ///< The number of successful steps since decreased precision.

			mutable mpfr endtime_highest_precision_;

		public:
			// functions offered for observers.
			template<typename RT>
			RT CondNum() const { return std::get<RT>(condition_number_estimate_);}

			template<typename RT>
			RT NormJ() const { return std::get<RT>(norm_J_);}

			template<typename RT>
			RT NormJInv() const { return std::get<RT>(norm_J_inverse_);}

			unsigned CurrentPrecision() const override
			{
				return current_precision_;
			}
		}; // re: class Tracker

	} // namespace tracking
} // namespace bertini


#endif



