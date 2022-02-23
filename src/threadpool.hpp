#include<thread>
#include<queue>
#include<vector>
#include<future>
#include<memory>
#include<functional>

using namespace std;

#define THREAD_POOL

class ThreadPool{
	struct TaskBlock{
		TaskBlock(): isTermination(true){}
		TaskBlock(const shared_ptr<const vector<function<score_t(const int)> > > &pFuncs, promise<vector<score_t> > results, future<TaskBlock> next):
			isTermination(false), pFuncs(pFuncs), results(move(results)), next(move(next)){}
		
		bool isTermination;
		const shared_ptr<const vector<function<score_t(const int)> > > pFuncs;
		promise<vector<score_t> > results;
		future<TaskBlock> next;
	};
	
	int nThreads;
	vector<promise<TaskBlock> > tasks;
	vector<function<score_t(const int)> > funcs;
	queue<score_t> sums;

	static void worker(future<TaskBlock> next, const int id){
		while (true){
			TaskBlock task = next.get();
			if (task.isTermination) break;
			vector<score_t> results;
			for (const function<score_t(const int)> &f: *(task.pFuncs)){
				results.push_back(f(id));
			}
			task.results.set_value(move(results));
			next = move(task.next);
		}
	}
	
	#ifdef USE_CUDA
	void work(){
		for (const function<score_t(const int)> &f: funcs) f(0);
		funcs.clear();
		gpuWork(sums);
	}
	#else
	void work(){
		const shared_ptr<const vector<function<score_t(const int)> > > pFuncs(new const vector<function<score_t(const int)> >(move(funcs)));
		funcs = vector<function<score_t(const int)> >();
		vector<future<vector<score_t> > > results;
		for (promise<TaskBlock> &task: tasks){
			promise<TaskBlock> curTask = move(task);
			task = promise<TaskBlock>();
			promise<vector<score_t> > res;
			results.push_back(move(res.get_future()));
			curTask.set_value(TaskBlock(pFuncs, move(res), task.get_future()));
		}
		vector<vector<score_t> > resultValues;
		resultValues.emplace_back();
		for (const function<score_t(const int)> &f: *(pFuncs)){
			resultValues[0].push_back(f(0));
		}
		for (future<vector<score_t> > &result: results){
			resultValues.push_back(result.get());
		}
		for (int i = 0; i < resultValues[0].size(); i++){
			score_t s = 0;
			for (vector<score_t> &rv: resultValues) s += rv[i];
			sums.push(s);
		}
	}
	#endif

public:
	#ifdef USE_CUDA
	static function<void(queue<score_t>&)> gpuWork;
	#endif

	ThreadPool(){}

	ThreadPool(int n): nThreads(n - 1), tasks(n - 1){
		for (int id = 1; id < n; id++){
			thread(worker, tasks[id - 1].get_future(), id).detach();
		}
	}

	void initialize(int n){
		nThreads = n - 1;
		tasks.resize(n - 1);
		for (int id = 1; id < n; id++){
			thread(worker, tasks[id - 1].get_future(), id).detach();
		}
	}
	
	~ThreadPool(){
		for (promise<TaskBlock> &task: tasks){
			task.set_value(TaskBlock());
		}
	}
	
	void push(function<score_t(const int)> f){
		funcs.push_back(f);
	}
	
	score_t pop(){
		if (sums.empty()) work();
		score_t v = sums.front();
		sums.pop();
		return v;
	}
};

#ifdef USE_CUDA
function<void(queue<score_t>&)> ThreadPool::gpuWork;
#endif