const buildForm = document.getElementById("stratifier-form");
const buildError = document.getElementById("build-error");
const buildButton = buildForm.querySelector("button[type=submit]");
let requestKey = null;
let submittedPrompt = "";

buildForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  const prompt = buildForm.elements.prompt.value.trim();
  if (!requestKey || submittedPrompt !== prompt) {
    requestKey = crypto.randomUUID();
    submittedPrompt = prompt;
  }
  buildButton.disabled = true;
  buildButton.textContent = "Starting build...";
  buildError.hidden = true;
  try {
    const response = await fetch(buildForm.action, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ prompt, request_key: requestKey }),
    });
    const result = await response.json();
    if (!response.ok) throw new Error(result.error || "The server could not start this build.");
    window.location.assign(`/stratifiers?job=${encodeURIComponent(result.id)}`);
  } catch (error) {
    buildError.textContent = error instanceof SyntaxError
      ? "The server is restarting or unavailable. Retry in a moment; duplicate submissions are prevented."
      : error.message;
    buildError.hidden = false;
    buildButton.disabled = false;
    buildButton.textContent = "Find dataset and build";
  }
});

document.querySelectorAll(".retry-build").forEach((button) => {
  button.addEventListener("click", () => {
    requestKey = null;
    buildForm.elements.prompt.value = button.dataset.prompt;
    buildForm.requestSubmit();
  });
});

async function pollJob(element) {
  try {
    const response = await fetch(`/api/stratifier-jobs/${element.dataset.jobId}`, { cache: "no-store" });
    if (!response.ok) throw new Error("Status temporarily unavailable");
    const job = await response.json();
    element.querySelector(".job-stage").textContent = job.stage;
    if (job.status === "complete" || job.status === "failed") {
      window.location.reload();
      return;
    }
  } catch (error) {
    element.querySelector(".job-stage").textContent = "Reconnecting to build status...";
  }
  window.setTimeout(() => pollJob(element), 2500);
}

document.querySelectorAll('.build-job[data-status="queued"], .build-job[data-status="running"]')
  .forEach(pollJob);
