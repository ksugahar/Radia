// core_stubs.cpp — Minimal stubs replacing taskmanager.cpp, profiler.cpp,
// paje_trace.cpp to avoid static initialization issues in Cubit DLL context.
// Single-threaded only.

#include <core/ngcore.hpp>

namespace ngcore
{

// ============================================================
// TaskManager — single-threaded stub (no thread pool)
// ============================================================
bool TaskManager::use_paje_trace = false;
int TaskManager::max_threads = 1;
int TaskManager::num_threads = 1;
thread_local int TaskManager::thread_id = 0;

const std::function<void(TaskInfo&)> * TaskManager::func = nullptr;
const std::function<void()> * TaskManager::startup_function = nullptr;
const std::function<void()> * TaskManager::cleanup_function = nullptr;

std::atomic<int> TaskManager::ntasks{0};
Exception * TaskManager::ex = nullptr;
std::atomic<int> TaskManager::jobnr{0};
std::atomic<int> TaskManager::complete[8] = {};
std::atomic<int> TaskManager::done{0};
std::atomic<int> TaskManager::active_workers{0};
std::atomic<int> TaskManager::workers_on_node[8] = {};
int TaskManager::sleep_usecs = 1000;
bool TaskManager::sleep = false;
TaskManager::NodeData *TaskManager::nodedata[8] = {};
int TaskManager::num_nodes = 1;

TaskManager::TaskManager() { num_threads = 1; }
TaskManager::~TaskManager() {}
void TaskManager::StartWorkers() {}
void TaskManager::StopWorkers() {}
// SuspendWorkers, ResumeWorkers, GetMaxThreads are inline in header
void TaskManager::SetNumThreads(int) { num_threads = 1; }
int TaskManager::GetThreadId() { return thread_id; }
void TaskManager::CreateJob(const std::function<void(TaskInfo&)> &f, int antasks) {
  // Single-threaded fallback: run all tasks sequentially
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
std::array<size_t, NgProfiler::SIZE> NgProfiler::dummy_thread_times = {};
size_t * NgProfiler::thread_times = NgProfiler::dummy_thread_times.data();
std::array<size_t, NgProfiler::SIZE> NgProfiler::dummy_thread_flops = {};
size_t * NgProfiler::thread_flops = NgProfiler::dummy_thread_flops.data();
std::shared_ptr<Logger> NgProfiler::logger;

NgProfiler::NgProfiler() {}
NgProfiler::~NgProfiler() {}
void NgProfiler::Print(FILE*) {}
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
