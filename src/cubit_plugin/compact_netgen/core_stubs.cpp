// core_stubs.cpp — Minimal stubs replacing taskmanager.cpp, profiler.cpp,
// paje_trace.cpp to avoid static initialization issues in Cubit DLL context.
// Single-threaded only.

#include <core/ngcore.hpp>

namespace ngcore
{

// ============================================================
// TaskManager — single-threaded stub (no thread pool)
// ============================================================
// Netgen >= 6.2.2607 keeps the task manager per instance behind a
// thread_local pointer.  No instance is ever created here, so every
// parallel region runs inline on the calling thread.
thread_local TaskManager * task_manager = nullptr;
TaskManager * GetTaskManager() { return task_manager; }

bool TaskManager::use_paje_trace = false;
int TaskManager::max_threads = 1;
thread_local int TaskManager::thread_id = 0;
// -1 disables timers: NgProfiler then never touches per-worker storage.
thread_local int TaskManager::timer_thread_id = -1;
thread_local WorkerData * TaskManager::worker_data = nullptr;

WorkerData::~WorkerData() {}
TaskManager::TaskManager(int) : taskqueue_ptr(nullptr), num_nodes(1), num_threads(1) {}
TaskManager::~TaskManager() {}
void TaskManager::StartWorkers() {}
void TaskManager::StopWorkers() {}
void TaskManager::SetNumThreads(int) {}
int TaskManager::GetThreadId() { return thread_id; }
int TaskManager::GetTimerThreadId() { return timer_thread_id; }
WorkerData * TaskManager::GetWorkerData() { return worker_data; }
bool TaskManager::ProcessTask() { return false; }
void TaskManager::Loop(int) {}
std::list<std::tuple<std::string, double>> TaskManager::Timing() { return {}; }
void RunWithTaskManager(std::function<void()> alg) { alg(); }
int EnterTaskManager(int) { return 0; }
void ExitTaskManager(int) {}

void TaskManager::CreateJob(const std::function<void(TaskInfo&)> &f, int antasks) {
  // Single-threaded fallback: run all tasks sequentially.  antasks = -1
  // means "one task per thread", i.e. one task here.
  if (antasks < 1) antasks = 1;
  for (int i = 0; i < antasks; i++) {
    TaskInfo ti;
    ti.task_nr = i;
    ti.ntasks = antasks;
    ti.thread_nr = 0;
    ti.nthreads = 1;
    f(ti);
  }
}

// ============================================================
// NgProfiler — no-op stub (no timers vector allocation)
// ============================================================
std::vector<NgProfiler::TimerVal> NgProfiler::timers;
std::string NgProfiler::filename;
std::shared_ptr<Logger> NgProfiler::logger;

NgProfiler::NgProfiler() {}
NgProfiler::~NgProfiler() {}
void NgProfiler::Print(FILE*) {}
void NgProfiler::Reset() {}
int NgProfiler::CreateTimer(const std::string &) {
  // Must return a valid index into timers vector
  timers.push_back(TimerVal());
  return (int)timers.size() - 1;
}

// ============================================================
// PajeTrace — no-op stubs
// ============================================================
PajeTrace *trace = nullptr;
std::vector<PajeTrace::MemoryEvent> PajeTrace::memory_events;
size_t PajeTrace::max_tracefile_size = 0;
bool PajeTrace::trace_thread_counter = false;
bool PajeTrace::trace_threads = false;
bool PajeTrace::mem_tracing_enabled = false;
bool PajeTrace::write_paje_file = false;
PajeTrace::PajeTrace(int, std::string) {}
PajeTrace::~PajeTrace() {}
void PajeTrace::StopTracing() {}
void PajeTrace::WritePajeFile(const std::string&) {}

// ============================================================
// Logging stubs (replace logging.cpp to avoid static init)
// ============================================================
// Must not be nullptr — Netgen code dereferences *testout unconditionally.
// Use std::clog (pre-existing global, no file-scope static init needed).
std::ostream* testout = &std::clog;
}  // close ngcore namespace

namespace netgen {
std::ostream* myerr = &std::cerr;
std::ostream* mycout = &std::cout;
void MyError(const char* ch) { if (myerr) (*myerr) << ch << std::endl; }
}

namespace ngcore {  // reopen ngcore
level::level_enum Logger::global_level = level::warn;
void Logger::log(level::level_enum level, std::string && s) {
  if (level >= global_level)
    std::clog << s << '\n';
}

} // namespace ngcore

// spdlog stub
namespace spdlog { class logger { public: logger() = default; }; }

namespace ngcore {
std::shared_ptr<Logger> GetLogger(const std::string&) {
  return std::make_shared<Logger>(std::make_shared<spdlog::logger>());
}
void SetLoggingLevel(level::level_enum level, const std::string&) {
  Logger::SetGlobalLoggingLevel(level);
}
void AddFileSink(const std::string&, level::level_enum, const std::string&) {}
void AddConsoleSink(level::level_enum, const std::string&) {}
void ClearLoggingSinks(const std::string&) {}
void FlushOnLoggingLevel(level::level_enum, const std::string&) {}
} // namespace ngcore
