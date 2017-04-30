//This file is part of Bertini 2.0.
//
// src/tracking/fixed_precision_tracker.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// src/tracking/fixed_precision_tracker.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  src/tracking/fixed_precision_tracker.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M




#include "bertini2/trackers/fixed_precision_tracker.hpp"

namespace bertini{
namespace tracking{

template<class DerivedT>
FixedPrecisionTracker<DerivedT>::FixedPrecisionTracker(System const& sys) : Base(sys){}


template<class DerivedT>
SuccessCode FixedPrecisionTracker<DerivedT>::TrackerIteration() const
{
	static_assert(std::is_same<	typename Eigen::NumTraits<RT>::Real, 
              				typename Eigen::NumTraits<CT>::Real>::value,
              				"underlying complex type and the type for comparisons must match");

	this->NotifyObservers(NewStep<EmitterType >(*this));

	Vec<CT>& predicted_space = std::get<Vec<CT> >(this->temporary_space_); // this will be populated in the Predict step
	Vec<CT>& current_space = std::get<Vec<CT> >(this->current_space_); // the thing we ultimately wish to update
	CT current_time = CT(this->current_time_);
	CT delta_t = CT(this->delta_t_);

	SuccessCode predictor_code = Predict(predicted_space, current_space, current_time, delta_t);

	if (predictor_code!=SuccessCode::Success)
	{
		this->NotifyObservers(FirstStepPredictorMatrixSolveFailure<EmitterType >(*this));

		this->next_stepsize_ = Get<Stepping>().step_size_fail_factor*this->current_stepsize_;

		UpdateStepsize();

		return predictor_code;
	}

	this->NotifyObservers(SuccessfulPredict<EmitterType , CT>(*this, predicted_space));

	Vec<CT>& tentative_next_space = std::get<Vec<CT> >(this->tentative_space_); // this will be populated in the Correct step

	CT tentative_next_time = current_time + delta_t;

	SuccessCode corrector_code = Correct(tentative_next_space,
										 predicted_space,
										 tentative_next_time);

	if (corrector_code == SuccessCode::GoingToInfinity)
	{
		// there is no corrective action possible...
		return corrector_code;
	}
	else if (corrector_code!=SuccessCode::Success)
	{
		this->NotifyObservers(CorrectorMatrixSolveFailure<EmitterType >(*this));

		this->next_stepsize_ = Get<Stepping>().step_size_fail_factor*this->current_stepsize_;
		UpdateStepsize();

		return corrector_code;
	}

	
	this->NotifyObservers(SuccessfulCorrect<EmitterType , CT>(*this, tentative_next_space));

	// copy the tentative vector into the current space vector;
	current_space = tentative_next_space;
	return SuccessCode::Success;
}

template class FixedPrecisionTracker<DoublePrecisionTracker>;
template class FixedPrecisionTracker<MultiplePrecisionTracker>;
//////////
//
// fixed double precision
//
//////////


DoublePrecisionTracker::DoublePrecisionTracker(class System const& sys) : FixedPrecisionTracker<DoublePrecisionTracker>(sys)
{	}

unsigned DoublePrecisionTracker::CurrentPrecision() const
{
	return DoublePrecision();
}

SuccessCode DoublePrecisionTracker::TrackerLoopInitialization(BaseComplexType const& start_time,
                               BaseComplexType const& end_time,
							   Vec<BaseComplexType> const& start_point) const
{
	this->NotifyObservers(Initializing<EmitterType,BaseComplexType>(*this,start_time, end_time, start_point));

	// set up the master current time and the current step size
	this->current_time_ = start_time;
	this->endtime_ = end_time;
	std::get<Vec<BaseComplexType> >(this->current_space_) = start_point;
	if (this->reinitialize_stepsize_)
		this->SetStepSize(min(Get<Stepping>().initial_step_size,abs(start_time-end_time)/Get<Stepping>().min_num_steps));

	ResetCounters();

	return SuccessCode::Success;
}







/////////////
//
//  fixed multiple precision
//
///////////////


MultiplePrecisionTracker::MultiplePrecisionTracker(class System const& sys) : FixedPrecisionTracker<MultiplePrecisionTracker>(sys), precision_(DefaultPrecision())
{	}

unsigned MultiplePrecisionTracker::CurrentPrecision() const
{
	return precision_;
}


SuccessCode MultiplePrecisionTracker::TrackerLoopInitialization(BaseComplexType const& start_time,
                               BaseComplexType const& end_time,
							   Vec<BaseComplexType> const& start_point) const
{

	if (start_point(0).precision()!=DefaultPrecision())
	{
		std::stringstream err_msg;
		err_msg << "start point for fixed multiple precision tracker has differing precision from default (" << start_point(0).precision() << "!=" << DefaultPrecision() << "), tracking cannot start";
		throw std::runtime_error(err_msg.str());
	}

	if (start_point(0).precision()!=CurrentPrecision())
	{
		std::stringstream err_msg;
		err_msg << "start point for fixed multiple precision tracker has differing precision from tracker's precision (" << start_point(0).precision() << "!=" << CurrentPrecision() << "), tracking cannot start";
		throw std::runtime_error(err_msg.str());
	}

	if (DefaultPrecision()!=CurrentPrecision())
	{
		std::stringstream err_msg;
		err_msg << "current default precision differs from tracker's precision (" << DefaultPrecision() << "!=" << CurrentPrecision() << "), tracking cannot start";
		throw std::runtime_error(err_msg.str());
	}


	this->NotifyObservers(Initializing<EmitterType,BaseComplexType>(*this,start_time, end_time, start_point));

	// set up the master current time and the current step size
	this->current_time_ = start_time;
	this->endtime_ = end_time;
	std::get<Vec<BaseComplexType> >(this->current_space_) = start_point;
	if (this->reinitialize_stepsize_)
		this->SetStepSize(min(Get<Stepping>().initial_step_size,mpfr_float(abs(start_time-end_time)/Get<Stepping>().min_num_steps)));

	ResetCounters();

	return SuccessCode::Success;
}

bool MultiplePrecisionTracker::PrecisionSanityCheck() const
{	
	return GetSystem().precision() == precision_ &&
			DefaultPrecision()==precision_ && 
			std::get<Vec<mpfr> >(current_space_)(0).precision() == precision_ &&
			std::get<Vec<mpfr> >(tentative_space_)(0).precision() == precision_ &&
			std::get<Vec<mpfr> >(temporary_space_)(0).precision() == precision_ &&
			std::get<mpfr_float>(norm_delta_z_).precision() == precision_ &&
			std::get<mpfr_float>(condition_number_estimate_).precision() == precision_ &&
			std::get<mpfr_float>(error_estimate_).precision() == precision_ &&
			std::get<mpfr_float>(norm_J_).precision() == precision_ &&
			std::get<mpfr_float>(norm_J_inverse_).precision() == precision_ &&
			std::get<mpfr_float>(size_proportion_).precision() == precision_ && 
			Precision(this->endtime_)==precision_
			        ;				
}


}// re: tracking
}// re: bertini





