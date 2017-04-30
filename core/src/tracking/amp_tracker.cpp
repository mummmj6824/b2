//This file is part of Bertini 2.0.
//
// src/tracking/amp_tracker.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// src/tracking/amp_tracker.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  src/tracking/amp_tracker.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M




#include "bertini2/trackers/amp_tracker.hpp"

namespace bertini{
namespace tracking{


void MinimizeTrackingCost(unsigned & new_precision, mpfr_float & new_stepsize, 
				  unsigned min_precision, mpfr_float const& old_stepsize,
				  unsigned max_precision, mpfr_float const& max_stepsize,
				  mpfr_float const& criterion_B_rhs,
				  unsigned num_newton_iterations,
				  unsigned predictor_order)
{
	mpfr_float min_cost = Eigen::NumTraits<mpfr_float>::highest();
	new_precision = MaxPrecisionAllowed()+1; // initialize to an impossible value.
	new_stepsize = old_stepsize; // initialize to original step size.

	auto minimizer_routine = 
		[&min_cost, &new_stepsize, &new_precision, criterion_B_rhs, num_newton_iterations, predictor_order, max_stepsize](unsigned candidate_precision)
		{
			mpfr_float candidate_stepsize = min(StepsizeSatisfyingCriterionB(candidate_precision, criterion_B_rhs, num_newton_iterations, predictor_order),
			                              max_stepsize);

			mpfr_float current_cost = ArithmeticCost(candidate_precision) / abs(candidate_stepsize);

			if (current_cost < min_cost)
			{
				min_cost = current_cost;
				new_stepsize = candidate_stepsize;
				new_precision = candidate_precision;
			}
		};

	unsigned lowest_mp_precision_to_test = min_precision;

	if (min_precision<=DoublePrecision())
		minimizer_routine(DoublePrecision());			


	if (lowest_mp_precision_to_test < LowestMultiplePrecision())
		lowest_mp_precision_to_test = LowestMultiplePrecision();
	else
		lowest_mp_precision_to_test = (lowest_mp_precision_to_test/PrecisionIncrement()) * PrecisionIncrement(); // use integer arithmetic to round.


	if (max_precision < LowestMultiplePrecision())
		max_precision = LowestMultiplePrecision();
	else
		max_precision = (max_precision/PrecisionIncrement()) * PrecisionIncrement(); // use integer arithmetic to round.

	for (unsigned candidate_precision = lowest_mp_precision_to_test; candidate_precision <= max_precision; candidate_precision+=PrecisionIncrement())
		minimizer_routine(candidate_precision);
}


SuccessCode AMPTracker::TrackerLoopInitialization(mpfr const& start_time,
                               mpfr const& end_time,
							   Vec<mpfr> const& start_point) const
{
	#ifndef BERTINI_DISABLE_ASSERTS
	assert(
	        (!preserve_precision_ 
	         || 
	         -log10(tracking_tolerance_) <= start_point(0).precision())
	         && "when tracking a path, either preservation of precision must be turned off (so precision can be higher at the end of tracking), or the initial precision must be high enough to support the resulting points to the desired tolerance"
	         );
	#endif

	NotifyObservers(Initializing<AMPTracker,mpfr>(*this,start_time, end_time, start_point));

	initial_precision_ = Precision(start_point(0));
	DefaultPrecision(initial_precision_);
	// set up the master current time and the current step size
	current_time_.precision(initial_precision_);
	current_time_ = start_time;

	endtime_highest_precision_.precision(initial_precision_);
	endtime_highest_precision_ = end_time;

	endtime_.precision(initial_precision_);
	endtime_ = end_time;

	current_stepsize_.precision(initial_precision_);
	if (reinitialize_stepsize_)
	{
		mpfr_float segment_length = abs(start_time-end_time)/Get<Stepping>().min_num_steps;
		SetStepSize(min(Get<Stepping>().initial_step_size,segment_length));
	}

	// populate the current space value with the start point, in appropriate precision
	if (initial_precision_==DoublePrecision())
		MultipleToDouble( start_point);
	else
		MultipleToMultiple(initial_precision_, start_point);
	
	ChangePrecision<upsample_refine_off>(initial_precision_);
	
	ResetCounters();

	ChangePrecision<upsample_refine_off>(start_point(0).precision());

	return InitialRefinement();
}



void AMPTracker::ResetCounters() const
{
	Tracker::ResetCountersBase();
	num_precision_decreases_ = 0;
	num_successful_steps_since_stepsize_increase_ = 0;
	num_successful_steps_since_precision_decrease_ = 0;
	// initialize to the frequency so guaranteed to compute it the first try 	
	num_steps_since_last_condition_number_computation_ = this->Get<Stepping>().frequency_of_CN_estimation;
}

SuccessCode AMPTracker::InitialRefinement() const
{
	SuccessCode initial_refinement_code = RefineStoredPoint();
	if (initial_refinement_code!=SuccessCode::Success)
	{
		do {
			if (current_precision_ > Get<PrecConf>().maximum_precision)
			{
				NotifyObservers(SingularStartPoint<EmitterType>(*this));
				return SuccessCode::SingularStartPoint;
			}

			if (current_precision_==DoublePrecision())
				initial_refinement_code = ChangePrecision<upsample_refine_on>(LowestMultiplePrecision());
			else
				initial_refinement_code = ChangePrecision<upsample_refine_on>(current_precision_+PrecisionIncrement());
		}
		while (initial_refinement_code!=SuccessCode::Success);
	}
	return SuccessCode::Success;
}	


SuccessCode AMPTracker::PreIterationCheck() const
{
	if (num_successful_steps_taken_ >= Get<Stepping>().max_num_steps)
		return SuccessCode::MaxNumStepsTaken;
	if (current_stepsize_ < Get<Stepping>().min_step_size)
		return SuccessCode::MinStepSizeReached;
	if (current_precision_ > Get<PrecConf>().maximum_precision)
		return SuccessCode::MaxPrecisionReached;

	return SuccessCode::Success;
}


void AMPTracker::PostTrackCleanup() const
{
	if (preserve_precision_)
		ChangePrecision(initial_precision_);
	NotifyObservers(TrackingEnded<EmitterType>(*this));
}


void AMPTracker::CopyFinalSolution(Vec<mpfr> & solution_at_endtime) const
{

	// the current precision is the precision of the output solution point.
	if (current_precision_==DoublePrecision())
	{
		unsigned num_vars = GetSystem().NumVariables();
		solution_at_endtime.resize(num_vars);
		for (unsigned ii=0; ii<num_vars; ii++)
			solution_at_endtime(ii) = mpfr(std::get<Vec<dbl> >(current_space_)(ii));
	}
	else
	{
		unsigned num_vars = GetSystem().NumVariables();
		solution_at_endtime.resize(num_vars);
		for (unsigned ii=0; ii<num_vars; ii++)
		{
			solution_at_endtime(ii).precision(current_precision_);
			solution_at_endtime(ii) = std::get<Vec<mpfr> >(current_space_)(ii);
		}
	}
}




SuccessCode AMPTracker::TrackerIteration() const
{
	if (current_precision_==DoublePrecision())
		return TrackerIteration<dbl, double>();
	else
		return TrackerIteration<mpfr, mpfr_float>();
}



template <typename ComplexType, typename RealType>
SuccessCode AMPTracker::TrackerIteration() const
{	
	static_assert(std::is_same<	typename Eigen::NumTraits<RealType>::Real, 
              				typename Eigen::NumTraits<ComplexType>::Real>::value,
              				"underlying complex type and the type for comparisons must match");

	#ifndef BERTINI_DISABLE_ASSERTS
	assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
	#endif

	NotifyObservers(NewStep<EmitterType>(*this));

	Vec<ComplexType>& predicted_space = std::get<Vec<ComplexType> >(temporary_space_); // this will be populated in the Predict step
	Vec<ComplexType>& current_space = std::get<Vec<ComplexType> >(current_space_); // the thing we ultimately wish to update
	ComplexType current_time = ComplexType(current_time_);
	ComplexType delta_t = ComplexType(delta_t_);

	SuccessCode predictor_code = Predict<ComplexType, RealType>(predicted_space, current_space, current_time, delta_t);
	if (predictor_code==SuccessCode::MatrixSolveFailureFirstPartOfPrediction)
	{
		NotifyObservers(FirstStepPredictorMatrixSolveFailure<EmitterType>(*this));
		next_stepsize_ = current_stepsize_;

		if (current_precision_==DoublePrecision())
			next_precision_ = LowestMultiplePrecision();
		else
			next_precision_ = current_precision_+(1+num_consecutive_failed_steps_) * PrecisionIncrement();

		UpdatePrecisionAndStepsize();

		return predictor_code;
	}
	else if (predictor_code==SuccessCode::MatrixSolveFailure)
	{
		NotifyObservers(PredictorMatrixSolveFailure<EmitterType>(*this));
		NewtonConvergenceError();// decrease stepsize, and adjust precision as necessary
		return predictor_code;
	}	
	else if (predictor_code==SuccessCode::HigherPrecisionNecessary)
	{	
		NotifyObservers(PredictorHigherPrecisionNecessary<EmitterType>(*this));
		AMPCriterionError<ComplexType, RealType>();
		return predictor_code;
	}


	NotifyObservers(SuccessfulPredict<AMPTracker, ComplexType>(*this, predicted_space));

	Vec<ComplexType>& tentative_next_space = std::get<Vec<ComplexType> >(tentative_space_); // this will be populated in the Correct step

	ComplexType tentative_next_time = current_time + delta_t;

	SuccessCode corrector_code = Correct<ComplexType, RealType>(tentative_next_space,
										 predicted_space,
										 tentative_next_time);

	if (corrector_code==SuccessCode::MatrixSolveFailure || corrector_code==SuccessCode::FailedToConverge)
	{
		NotifyObservers(CorrectorMatrixSolveFailure<EmitterType>(*this));
		NewtonConvergenceError();
		return corrector_code;
	}
	else if (corrector_code == SuccessCode::HigherPrecisionNecessary)
	{
		NotifyObservers(CorrectorHigherPrecisionNecessary<EmitterType>(*this));
		AMPCriterionError<ComplexType, RealType>();
		return corrector_code;
	}
	else if (corrector_code == SuccessCode::GoingToInfinity)
	{
		// there is no corrective action possible...
		return corrector_code;
	}

	NotifyObservers(SuccessfulCorrect<AMPTracker, ComplexType>(*this, tentative_next_space));

	// copy the tentative vector into the current space vector;
	current_space = tentative_next_space;
	return AdjustAMPStepSuccess<ComplexType,RealType>();
}


template
SuccessCode AMPTracker::TrackerIteration<mpfr, mpfr_float>() const;

template
SuccessCode AMPTracker::TrackerIteration<dbl, double>() const;





template <typename ComplexType, typename RealType>
SuccessCode AMPTracker::AdjustAMPStepSuccess() const
{
	mpfr_float min_stepsize = current_stepsize_ * Get<Stepping>().step_size_fail_factor;
	mpfr_float max_stepsize = min(mpfr_float(current_stepsize_ * Get<Stepping>().step_size_success_factor), Get<Stepping>().max_step_size);


	unsigned min_precision = MinRequiredPrecision_BCTol<ComplexType, RealType>();
	unsigned max_precision = max(min_precision,current_precision_);

	if (num_successful_steps_since_stepsize_increase_ < Get<Stepping>().consecutive_successful_steps_before_stepsize_increase)
		max_stepsize = current_stepsize_; // disallow stepsize changing 


	if ( (num_successful_steps_since_precision_decrease_ < Get<PrecConf>().consecutive_successful_steps_before_precision_decrease)
	    ||
	    (num_precision_decreases_ >= Get<PrecConf>().max_num_precision_decreases))
		min_precision = max(min_precision, current_precision_); // disallow precision changing 


	MinimizeTrackingCost(next_precision_, next_stepsize_, 
				min_precision, min_stepsize,
				max_precision, max_stepsize,
				B_RHS<ComplexType, RealType>(),
				Get<Newton>().max_num_newton_iterations,
				predictor_order_);


	if ( (next_stepsize_ > current_stepsize_) || (next_precision_ < current_precision_) )
		num_successful_steps_since_stepsize_increase_ = 0;
	else
		num_successful_steps_since_stepsize_increase_++;
	

	if (next_precision_ < current_precision_)
	{ 
		num_precision_decreases_++;
		num_successful_steps_since_precision_decrease_ = 0;
	}
	else
		num_successful_steps_since_precision_decrease_ ++;

	return UpdatePrecisionAndStepsize();
}


template
SuccessCode AMPTracker::AdjustAMPStepSuccess<mpfr,mpfr_float>() const;

template
SuccessCode AMPTracker::AdjustAMPStepSuccess<dbl,double>() const;






void AMPTracker::NewtonConvergenceError() const
{
	next_precision_ = current_precision_;
	next_stepsize_ = Get<Stepping>().step_size_fail_factor*current_stepsize_;

	while (next_stepsize_ < MinStepSizeForPrecision(next_precision_))
	{
		if (next_precision_==DoublePrecision())
			next_precision_=LowestMultiplePrecision();
		else
			next_precision_+=PrecisionIncrement();
	}

	UpdatePrecisionAndStepsize();
}









template<typename ComplexType, typename RealType>
void AMPTracker::AMPCriterionError() const
{	
	unsigned min_next_precision;
	if (current_precision_==DoublePrecision())
		min_next_precision = LowestMultiplePrecision(); // precision increases
	else
		min_next_precision = current_precision_ + (1+num_consecutive_failed_steps_)*PrecisionIncrement(); // precision increases


	mpfr_float min_stepsize = MinStepSizeForPrecision(current_precision_);
	mpfr_float max_stepsize = current_stepsize_ * Get<Stepping>().step_size_fail_factor;  // Stepsize decreases.

	if (min_stepsize > max_stepsize)
	{
		// stepsizes are incompatible, must increase precision
		next_precision_ = min_next_precision;
		// decrease stepsize somewhat less than the fail factor
		next_stepsize_ = current_stepsize_ * (1+Get<Stepping>().step_size_fail_factor)/2;
	}
	else
	{
		mpfr_float criterion_B_rhs = B_RHS<ComplexType,RealType>();

		unsigned min_precision = max(min_next_precision,
		                             criterion_B_rhs.convert_to<unsigned int>(),
		                             DigitsC<ComplexType,RealType>(),
		                             MinDigitsForStepsizeInterval(min_stepsize, max_stepsize, current_time_),
		                             digits_final_
		                             );
		
		mpfr_float a(ceil(criterion_B_rhs - (predictor_order_+1)* -log10(max_stepsize)/Get<Newton>().max_num_newton_iterations));

		unsigned max_precision = max(min_precision, 
		                             a.convert_to<unsigned int>()
		                             );

		MinimizeTrackingCost(next_precision_, next_stepsize_, 
				min_precision, min_stepsize,
				Get<PrecConf>().maximum_precision, max_stepsize,
				criterion_B_rhs,
				Get<Newton>().max_num_newton_iterations,
				predictor_order_);
	}

	UpdatePrecisionAndStepsize();
}




template
void AMPTracker::AMPCriterionError<mpfr, mpfr_float>() const;

template
void AMPTracker::AMPCriterionError<dbl, double>() const;











void AMPTracker::DoubleToDouble(Vec<dbl> const& source_point) const
{	
	#ifndef BERTINI_DISABLE_ASSERTS
	assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
	#endif

	current_precision_ = DoublePrecision();
	DefaultPrecision(DoublePrecision());

	GetSystem().precision(16);

	std::get<Vec<dbl> >(current_space_) = source_point;
}






void AMPTracker::MultipleToDouble(Vec<mpfr> const& source_point) const
{	
	#ifndef BERTINI_DISABLE_ASSERTS
	assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
	#endif
	previous_precision_ = current_precision_;
	current_precision_ = DoublePrecision();
	DefaultPrecision(DoublePrecision());

	GetSystem().precision(DoublePrecision());

	if (std::get<Vec<dbl> >(current_space_).size()!=source_point.size())
		std::get<Vec<dbl> >(current_space_).resize(source_point.size());

	for (unsigned ii=0; ii<source_point.size(); ii++)
		std::get<Vec<dbl> >(current_space_)(ii) = dbl(source_point(ii));

	endtime_.precision(DoublePrecision());
}




void AMPTracker::DoubleToMultiple(unsigned new_precision, Vec<dbl> const& source_point) const
{	
	#ifndef BERTINI_DISABLE_ASSERTS
	assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
	assert(new_precision > DoublePrecision() && "must convert to precision higher than DoublePrecision when converting to multiple precision");
	#endif
	previous_precision_ = current_precision_;
	current_precision_ = new_precision;
	DefaultPrecision(new_precision);
	GetSystem().precision(new_precision);
	predictor_->ChangePrecision(new_precision);
	corrector_->ChangePrecision(new_precision);

	endtime_ = endtime_highest_precision_;
	endtime_.precision(new_precision);

	current_time_.precision(new_precision);

	if (std::get<Vec<mpfr> >(current_space_).size()!=source_point.size())
		std::get<Vec<mpfr> >(current_space_).resize(source_point.size());

	for (unsigned ii=0; ii<source_point.size(); ii++)
		std::get<Vec<mpfr> >(current_space_)(ii) = mpfr(source_point(ii));

	AdjustTemporariesPrecision(new_precision);

	#ifndef BERTINI_DISABLE_ASSERTS
	assert(std::get<Vec<mpfr> >(current_space_)(0).precision() == current_precision_ && "precision of time in mpfr doesn't match tracker");
	#endif
}




void AMPTracker::MultipleToMultiple(unsigned new_precision, Vec<mpfr> const& source_point) const
{	
	#ifndef BERTINI_DISABLE_ASSERTS
	assert(source_point.size() == GetSystem().NumVariables() && "source point for converting to multiple precision is not the same size as the number of variables in the system being solved.");
	assert(new_precision > DoublePrecision() && "must convert to precision higher than DoublePrecision when converting to multiple precision");
	#endif
	previous_precision_ = current_precision_;
	current_precision_ = new_precision;
	DefaultPrecision(new_precision);
	GetSystem().precision(new_precision);
	predictor_->ChangePrecision(new_precision);
	corrector_->ChangePrecision(new_precision);

	endtime_ = endtime_highest_precision_;
	endtime_.precision(new_precision);

	current_time_.precision(new_precision);

	if (std::get<Vec<mpfr> >(current_space_).size()!=source_point.size())
		std::get<Vec<mpfr> >(current_space_).resize(source_point.size());

	for (unsigned ii=0; ii<source_point.size(); ii++)
		std::get<Vec<mpfr> >(current_space_)(ii) = mpfr(source_point(ii));

	AdjustTemporariesPrecision(new_precision);

	#ifndef BERTINI_DISABLE_ASSERTS
	assert(std::get<Vec<mpfr> >(current_space_)(0).precision() == current_precision_ && "precision of time in mpfr doesn't match tracker");
	#endif
}




void AMPTracker::AdjustTemporariesPrecision(unsigned new_precision) const
{
	unsigned num_vars = GetSystem().NumVariables();

	//  the current_space value is adjusted in the appropriate ChangePrecision function
	std::get<Vec<mpfr> >(tentative_space_).resize(num_vars);
	for (unsigned ii = 0; ii < num_vars; ++ii)
		std::get<Vec<mpfr> >(tentative_space_)(ii).precision(new_precision);

	std::get<Vec<mpfr> >(temporary_space_).resize(num_vars);
	for (unsigned ii = 0; ii < num_vars; ++ii)
		std::get<Vec<mpfr> >(temporary_space_)(ii).precision(new_precision);

	std::get<mpfr_float>(condition_number_estimate_).precision(new_precision);
	std::get<mpfr_float>(error_estimate_).precision(new_precision);
	std::get<mpfr_float>(norm_J_).precision(new_precision);
	std::get<mpfr_float>(norm_J_inverse_).precision(new_precision);
	std::get<mpfr_float>(norm_delta_z_).precision(new_precision);
	std::get<mpfr_float>(size_proportion_).precision(new_precision);
}


bool AMPTracker::PrecisionSanityCheck() const
{	
	if (current_precision_==DoublePrecision())
	{
		return true;
	}
	else
	{
		assert(DefaultPrecision()==current_precision_ && "current precision differs from the default precision");

		return GetSystem().precision() == current_precision_ &&
				predictor_->precision() == current_precision_ &&
				std::get<Vec<mpfr> >(current_space_)(0).precision() == current_precision_ &&
				std::get<Vec<mpfr> >(tentative_space_)(0).precision() == current_precision_ &&
				std::get<Vec<mpfr> >(temporary_space_)(0).precision() == current_precision_ &&
				std::get<mpfr_float>(norm_delta_z_).precision() == current_precision_ &&
				std::get<mpfr_float>(condition_number_estimate_).precision() == current_precision_ &&
				std::get<mpfr_float>(error_estimate_).precision() == current_precision_ &&
				std::get<mpfr_float>(norm_J_).precision() == current_precision_ &&
				std::get<mpfr_float>(norm_J_inverse_).precision() == current_precision_ &&
				std::get<mpfr_float>(size_proportion_).precision() == current_precision_ &&
				Precision(endtime_) == current_precision_ && 
				Precision(current_time_) == current_precision_
				        ;
	}
	
}


}// re: tracking
}// re: bertini







