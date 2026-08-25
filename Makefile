VENV_BIN := .venv/bin

.DEFAULT_GOAL := help
.PHONY: help sync format format-check lint typecheck deps imports test docs docs-serve build carpets check clean

help:  ## Show developer commands
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

sync:  ## Synchronize the locked development environment
	uv sync --frozen --group dev --group docs

format:  ## Format and autofix source, tests, and benchmarks
	$(VENV_BIN)/ruff format src tests benchmarks
	$(VENV_BIN)/ruff check --fix src tests benchmarks

format-check:  ## Check formatting without changing files
	$(VENV_BIN)/ruff format --check src tests benchmarks

lint:  ## Run Ruff lint checks
	$(VENV_BIN)/ruff check src tests benchmarks

typecheck:  ## Run standard Pyright in strict mode
	$(VENV_BIN)/pyright

deps:  ## Validate dependency declarations
	$(VENV_BIN)/deptry .

imports:  ## Enforce directed package boundaries
	$(VENV_BIN)/lint-imports

test:  ## Run tests with branch coverage
	$(VENV_BIN)/pytest --cov --cov-branch

docs:  ## Build documentation and fail on warnings
	$(VENV_BIN)/mkdocs build --strict

docs-serve:  ## Serve documentation with live reload
	$(VENV_BIN)/mkdocs serve

build:  ## Build and validate source and wheel distributions
	uv build
	$(VENV_BIN)/twine check dist/*

carpets:  ## Report carpet diagnostics (never blocks; flags for triage, not verdicts)
	$(VENV_BIN)/carpet-scan src --tests tests --html build/carpet-report.html

check:  ## Run every merge-blocking quality gate
	uv lock --check
	$(MAKE) format-check lint typecheck deps imports test docs build

clean:  ## Remove generated build and quality artifacts
	$(VENV_BIN)/python -c "import shutil; [shutil.rmtree(path, ignore_errors=True) for path in ('build', 'dist', 'site', '.pytest_cache', '.ruff_cache')]"
