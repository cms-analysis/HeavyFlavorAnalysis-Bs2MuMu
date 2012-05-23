#include <RooStats/HypoTestInverterResult.h>

const static double sm=3.2e-9;

class MyResult : public RooStats::HypoTestInverterResult {
	public:
		MyResult(RooStats::HypoTestInverterResult *cpy) {
			Add(*cpy); // add this result
		}
		
		void rescale(double scale) {
			for (unsigned j = 0; j < fXValues.size(); j++) fXValues[j] *= scale;
		}
	private:
		ClassDef(MyResult,1)
};

RooStats::HypoTestInverterResult *convert(RooStats::HypoTestInverterResult *inres)
{
	MyResult *result = new MyResult(inres);
	result->rescale(sm);
	
	return result;
} // convert()
