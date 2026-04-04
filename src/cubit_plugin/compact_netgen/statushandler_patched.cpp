// statushandler_patched.cpp — Compact Netgen
// All file-scope statics converted to function-local statics (lazy init)
// to avoid DLL DllMain static initialization issues.

#include "array.hpp"
#include "statushandler.hpp"

namespace ngcore
{
  volatile multithreadt multithread;

  multithreadt :: multithreadt()
  {
    pause = 0;
    testmode = 0;
    redraw = 0;
    drawing = 0;
    terminate = 0;
    running = 0;
    percent = 0;
    task = "";
  }

  // Lazy-init wrappers (avoid file-scope static init in DLL context)
  static Array<std::string>& get_msgstatus_stack() {
    static Array<std::string> s(0);
    return s;
  }
  static Array<double>& get_threadpercent_stack() {
    static Array<double> s(0);
    return s;
  }
  static std::string& get_msgstatus() {
    static std::string s;
    return s;
  }

  void ResetStatus()
  {
    SetStatMsg("idle");
    get_msgstatus_stack().SetSize(0);
    get_threadpercent_stack().SetSize(0);
    multithread.percent = 100.;
  }

  void PushStatus(const std::string& s)
  {
    get_msgstatus_stack().Append(s);
    SetStatMsg(s);
    get_threadpercent_stack().Append(0);
  }

  void PopStatus()
  {
    auto& msgstatus_stack = get_msgstatus_stack();
    auto& threadpercent_stack = get_threadpercent_stack();
    if (msgstatus_stack.Size()) {
      if (msgstatus_stack.Size() > 1)
        SetStatMsg(msgstatus_stack[msgstatus_stack.Size()-2]);
      else
        SetStatMsg("");
      msgstatus_stack.DeleteLast();
      threadpercent_stack.DeleteLast();
      if (threadpercent_stack.Size() > 0)
        multithread.percent = threadpercent_stack.Last();
      else
        multithread.percent = 100.;
    }
  }

  void SetStatMsg(const std::string& s)
  {
    get_msgstatus() = s;
    multithread.task = get_msgstatus().c_str();
  }

  void SetThreadPercent(double percent)
  {
    multithread.percent = percent;
    if (get_threadpercent_stack().Size() > 0)
      get_threadpercent_stack().Last() = percent;
  }

  void GetStatus(std::string & s, double & percentage)
  {
    if (get_threadpercent_stack().Size() > 0)
      percentage = get_threadpercent_stack().Last();
    else
      percentage = multithread.percent;

    if (get_msgstatus_stack().Size())
      s = get_msgstatus_stack().Last();
    else
      s = "idle";
  }
}
